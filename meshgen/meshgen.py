from itertools import izip, repeat
import multiprocessing as mp
from operator import methodcaller
from os.path import abspath, join
import os, sys

import numpy as np

from GoTools import Point, WriteG2
from GoTools.SurfaceFactory import LoftCurves
import GeoUtils.Interpolate as ip
from GeoUtils.IO import Numberer, InputFile

from utils import *
from airfoil import AirFoil
from slice import Slice
from volumetricslice import VolumetricSlice


def fill_runner(af):
    """Produces slices from airfoils. Used for parallell computation."""
    return Slice.from_airfoil(af)


def progress(title, i, tot):
    """Progress output to stdout."""
    args = (title, i, tot, 100 * float(i) / tot)
    s = '\r%s: %i/%i (%.2f%%)' % args
    if i == tot:
        s += '\n'
    sys.stdout.write(s)
    sys.stdout.flush()


class MeshGen(object):
    """This is a utility class used for mesh generation. It will not do everything on its
    own. Methods must be called in the right order. Calling code can rely on just checking mesh_mode
    and length_mode."""

    def __init__(self, params):
        self.p = params

    def load_airfoils(self):
        """Loads airfoils from the wing definition file."""
        self.airfoils = [AirFoil.from_wingdef(self.p, wd) for wd in self.p.wingdef]
        self.p.dump_g2files('airfoils_raw', self.airfoils)

    def resample_airfoils(self):
        """Performs angular resampling of all airfoils."""
        for af in self.airfoils:
            af.resample()
        self.p.dump_g2files('airfoils_resampled_angular', self.airfoils)

    def resolve_join(self):
        """Resolves the join according to the join_adds and join_index parameters."""
        n = self.p.join_adds
        idx = self.p.join_index

        if n == 0:
            return

        # Interpolation upwards
        for _ in xrange(n):
            new = AirFoil.from_mean(self.airfoils[idx], self.airfoils[idx+1])
            self.airfoils.insert(idx + 1, new)

        # Interpolation downwards
        for i in xrange(n):
            new = AirFoil.from_mean(self.airfoils[idx+i-1], self.airfoils[idx+i])
            self.airfoils.insert(idx + i, new)

        # Remove the join
        del self.airfoils[idx + n]

        self.p.dump_g2files('airfoils_joined', self.airfoils)

    def resample_length(self):
        """Resamples airfoils in the length direction."""
        # Loft the curves to produce a wing surface
        curves = [af.curve for af in self.airfoils]
        zvals = [af.z() for af in self.airfoils]
        wing = LoftCurves(curves, zvals, order=4)

        # Interpolating curve for theta as a function of z
        thetas = ip.CubicCurve(x=[af.theta for af in self.airfoils],
                               t=zvals, boundary=ip.NATURAL)
        theta = lambda z: thetas.Evaluate(z)[0]

        # Get the resampled z-values
        resampler = getattr(self, '_resample_length_' + self.p.length_mode)
        new_zvals = resampler(zvals[0], zvals[-1])

        # Produce new airfoils from point-evaluating the wing
        self.airfoils = []
        kus, _ = wing.GetKnots()
        for z in new_zvals:
            pts = [wing.Evaluate(ku, z) for ku in kus]
            self.airfoils.append(AirFoil.from_pts(self.p, theta(z), pts))

        self.p.dump_g2files('airfoils_resampled_length', self.airfoils)

    def make_slices(self):
        """Turns airfoils into slices."""
        progress('Making slices', 0, len(self.airfoils))
        pool = mp.Pool(self.p.nprocs_mg)
        result = []
        for i, af in enumerate(pool.imap(fill_runner, self.airfoils)):
            result.append(af)
            progress('Making slices', i+1, len(self.airfoils))

        del self.airfoils
        self.slices = result
        self.p.dump_g2files('slices_raw', self.slices)

    def extend(self):
        """Extend all the slices according to parameters (behind, ahead, sides)."""
        for s in self.slices:
            s.extend()
        self.p.dump_g2files('slices_sides', self.slices)

    def extrude(self):
        """Performs extrusion of a single airfoil or slice."""
        zvals = np.linspace(0, self.p.length, self.p.n_length + 1)
        attr = 'slices' if hasattr(self, 'slices') else 'airfoils'
        root = getattr(self, attr)[0]
        setattr(self, attr, [root.translate(Point(0,0,z)) for z in zvals])

    def loft_slices(self):
        """Lofts slices together to produce volumetric slices."""
        self.zvals = [s.z() for s in self.slices]

        temp = VolumetricSlice.from_slices(self.slices)
        self.slices = temp.subdivide()
        self.p.dump_g2files('slices_volumetric', self.slices)

    def loft_blade(self):
        """Lofts airfoils together to produce a blade."""
        params = extend_knots(map(float, range(len(self.airfoils))))
        airfoils = list(self.airfoils)
        airfoils.insert(1, AirFoil.from_mean(airfoils[0], airfoils[1]))
        airfoils.insert(-1, AirFoil.from_mean(airfoils[-2], airfoils[-1]))

        self.blade = LoftCurves([af.curve for af in airfoils], params, 4)

    def subdivide_slices(self):
        """Subdivices all the slices."""
        for s in self.slices:
            s.subdivide()
        self.p.dump_g2files('slices_subdivided', self.slices)

    def lower_order(self):
        """Lowers order on the slices and/or the blade."""
        if hasattr(self, 'slices'):
            for s in self.slices:
                s.lower_order()

        if hasattr(self, 'blade'):
            lower_order(self.blade, self.p.order)

    def output(self):
        """Produces the final output."""
        if self.p.format == 'OpenFOAM' and self.p.mesh_mode in {'2d', 'semi3d'}:
            self.slices = [VolumetricSlice.from_slice(s) for s in self.slices]

        getattr(self, '_output_' + self.p.mesh_mode)()
        self.p.out_yaml()

    def _output_beam(self):
        """Writes a beam.g2 file."""
        if not hasattr(self, 'zvals'):
            self.zvals = [s.z() for s in self.slices]

        beam = ip.LinearCurve(pts=[Point(0,0,z) for z in self.zvals])
        WriteG2(join(self.p.out, 'beam.g2'), beam)

    def _output_blade(self):
        """Produces the final output of the blade."""
        WriteG2(self.p.out_path() + '.g2', self.blade)

    def _output_2d(self):
        """Produces the final output in 2D mode."""
        self.slices[0].output(self.p.out_path())

    def _output_semi3d(self):
        """Produces the final output in semi3D mode."""
        progress('Writing planes', 0, len(self.slices))
        for i, s in enumerate(self.slices):
            custom = 'slice-%03i' % (i+1)
            path = self.p.out_path(custom)

            # In OpenFOAM mode, we have to output to separate subfolders
            if self.p.format == 'OpenFOAM':
                try:
                    os.makedirs(path)
                except OSError:
                    pass
                path = join(path, 'out')

            s.output(path, custom)
            progress('Writing planes', i+1, len(self.slices))

        if self.p.format == 'IFEM':
            self._output_beam()

    def _output_3d(self):
        """Produces the final output in 3D mode."""
        path = self.p.out_path()
        n = Numberer()

        # Add all the patches and boundaries
        for s in self.slices:
            s.push_patches(n)
            s.push_topsets(n)
        self.slices[0].push_hub(n)
        self.slices[-1].push_antihub(n)

        if self.p.debug:
            n.WriteTopologySets(path)

        if self.p.walldistance:
            n.AddWallGroup('wing')

        # Renumber and final output
        n.Renumber(self.p.nprocs)
        n.WriteEverything(path, display=True)

        if self.p.format == 'OpenFOAM':
            convert_openfoam(path)
        if self.p.format == 'IFEM':
            self._output_beam()
            self.p.postprocess_xinp()

    def _resample_length_uniform(self, za, zb):
        """Uniform length resampling."""
        return np.linspace(za, zb, self.p.n_length + 1)

    def _resample_length_double(self, za, zb):
        """Double-sided length resampling."""
        n = self.p.n_length / 2
        dj, dt = self.p.d_join, self.p.d_tip

        r1, r2 = grading_double(zb - za, dj, dt, n, n)

        zvals = gradspace(za, dj, r1, n+1)
        zvals += gradspace(zb, -dt, r2, n)[::-1]
        return zvals

    def _resample_length_triple(self, za, zb):
        """Triple-sided length resampling."""
        join = self.p.wingdef[self.p.join_index].z
        dj, dt = self.p.d_join, self.p.d_tip
        nb, nl = self.p.n_base, (self.p.n_length - self.p.n_base) / 2

        rz = grading(join - za, dj, nb)
        zvals = gradspace(join, -dj, rz, nb + 1)[::-1]

        rz1, rz2 = grading_double(zb - zvals[-1], dj, dt, nl, nl)
        zvals += gradspace(zvals[-1], dj, rz1, nl+1)[1:]
        zvals += gradspace(zb, -dt, rz2, nl)[::-1]

        return zvals
