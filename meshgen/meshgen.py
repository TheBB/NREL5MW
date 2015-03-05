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
    return Slice(af)


def progress(title, i, tot):
    args = (title, i, tot, 100 * float(i) / tot)
    s = '\r%s: %i/%i (%.2f%%)' % args
    if i == tot:
        s += '\n'
    sys.stdout.write(s)
    sys.stdout.flush()


class MeshGen(object):

    def __init__(self, params):
        self.p = params


    def load_airfoils(self):
        self.airfoils = [AirFoil.from_wingdef(self.p, wd) for wd in self.p.wingdef]
        self.p.dump_g2files('airfoils_raw', self.airfoils)


    def resample_airfoils(self):
        for af in self.airfoils:
            af.resample()
        self.p.dump_g2files('airfoils_resampled_radial', self.airfoils)


    def resolve_join(self):
        n = self.p.join_adds
        idx = self.p.join_index

        if n == 0:
            return

        for _ in xrange(n):
            new = AirFoil.from_mean(self.airfoils[idx], self.airfoils[idx+1])
            self.airfoils.insert(idx + 1, new)

        for i in xrange(n):
            new = AirFoil.from_mean(self.airfoils[idx+i-1], self.airfoils[idx+i])
            self.airfoils.insert(idx + i, new)

        del self.airfoils[idx + n]
        self.p.dump_g2files('airfoils_joined', self.airfoils)


    def resample_length(self):
        curves = [af.curve for af in self.airfoils]
        zvals = [af.z() for af in self.airfoils]
        wing = LoftCurves(curves, zvals, order=4)
        thetas = ip.CubicCurve(x=[af.theta for af in self.airfoils],
                               t=zvals, boundary=ip.NATURAL)
        theta = lambda z: thetas.Evaluate(z)[0]

        if self.p.length_mode == 'uniform':
            zvals = self._resample_length_uniform(zvals[0], zvals[-1])
        elif self.p.length_mode == 'double':
            zvals = self._resample_length_double(zvals[0], zvals[-1])
        elif self.p.length_mode == 'triple':
            zvals = self._resample_length_triple(zvals[0], zvals[-1])

        self.airfoils = []
        kus, _ = wing.GetKnots()
        for z in zvals:
            pts = [wing.Evaluate(ku, z) for ku in kus]
            self.airfoils.append(AirFoil.from_pts(self.p, theta(z), pts))

        self.p.dump_g2files('airfoils_resampled_length', self.airfoils)


    def make_slices(self):
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
        for s in self.slices:
            s.extend()
        self.p.dump_g2files('slices_sides', self.slices)


    def extrude(self):
        zvals = np.linspace(0, self.p.length, self.p.n_length + 1)
        attr = 'slices' if hasattr(self, 'slices') else 'airfoils'
        root = getattr(self, attr)[0]
        setattr(self, attr, [root.translate(Point(0,0,z)) for z in zvals])


    def loft_slices(self):
        temp = VolumetricSlice.from_slices(self.slices)
        self.slices = temp.subdivide()
        self.p.dump_g2files('slices_volumetric', self.slices)


    def loft_blade(self):
        params = map(float, range(len(self.airfoils)))
        airfoils = list(self.airfoils)

        airfoils.insert(1, AirFoil.from_mean(airfoils[0], airfoils[1]))
        params.insert(1, 0.5)

        airfoils.insert(-1, AirFoil.from_mean(airfoils[-2], airfoils[-1]))
        params.insert(-1, (params[-2] + params[-1])/2)

        self.blade = LoftCurves([af.curve for af in airfoils], params, 4)


    def subdivide_slices(self):
        for s in self.slices:
            s.subdivide()
        self.p.dump_g2files('slices_subdivided', self.slices)


    def lower_order(self):
        if hasattr(self, 'slices'):
            for s in self.slices:
                s.lower_order()

        if hasattr(self, 'blade'):
            lower_order(self.blade, self.p.order)


    def output(self):
        if self.p.format == 'OpenFOAM' and self.p.mesh_mode in {'2d', 'semi3d'}:
            self.slices = [VolumetricSlice.from_slice(s) for s in self.slices]

        getattr(self, '_output_' + self.p.mesh_mode)()
        self.p.out_yaml()


    def _output_blade(self):
        path = abspath(join(self.p.out, self.p.out)) + '.g2'
        WriteG2(path, self.blade)


    def _output_2d(self):
        path = abspath(join(self.p.out, self.p.out))
        self.slices[0].output(path)


    def _output_semi3d(self):
        progress('Writing planes', 0, len(self.slices))
        for i, s in enumerate(self.slices):
            path = abspath(join(self.p.out, 'slice-%03i' % (i+1)))
            if self.p.format == 'OpenFOAM':
                try:
                    os.makedirs(path)
                except OSError:
                    pass
                path = join(path, 'out')
            s.output(path)
            progress('Writing planes', i+1, len(self.slices))

        zvals = [s.z() for s in self.slices]
        beam = ip.LinearCurve(pts=[Point(z,0,0) for z in zvals])
        WriteG2(join(self.p.out, 'beam.g2'), beam)


    def _output_3d(self):
        path = abspath(join(self.p.out, self.p.out))
        n = Numberer()

        for s in self.slices:
            s.push_patches(n)
            s.push_boundaries(n)
        self.slices[0].push_hub(n)
        self.slices[-1].push_antihub(n)

        if self.p.debug:
            n.WriteBoundaries(path)

        if self.p.walldistance:
            n.AddWallGroup('wing')

        n.Renumber(self.p.nprocs)
        n.WriteEverything(path, display=True)

        if self.p.format == 'OpenFOAM':
            convert_openfoam(path)


    def _resample_length_uniform(self, za, zb):
        return np.linspace(za, zb, self.p.n_length + 1)


    def _resample_length_double(self, za, zb):
        n = self.p.n_length / 2
        dj, dt = self.p.d_join, self.p.d_tip

        r1, r2 = grading_double(zb - za, dj, dt, n-1, n-1)

        zvals = gradspace(za, dj, 1./r1, n+1)
        zvals += gradspace(zb, -dt, 1./r2, n)[::-1]
        return zvals


    def _resample_length_triple(self, za, zb):
        join = self.p.wingdef[self.p.join_index].z
        dj, dt = self.p.d_join, self.p.d_tip
        nb, nl = self.p.n_base, (self.p.n_length - self.p.n_base) / 2

        rz = grading(join - za, dj, nb)
        zvals = gradspace(join, -dj, rz, nb + 1)[::-1]

        rz1, rz2 = grading_double(zb - zvals[-1], dj, dt, nl-1, nl-1)
        zvals += gradspace(zvals[-1], dj, 1./rz1, nl+1)[1:]
        zvals += gradspace(zb, -dt, 1./rz2, nl)[::-1]

        return zvals
