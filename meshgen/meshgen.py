from itertools import izip, repeat
import multiprocessing as mp
from operator import methodcaller
from os.path import abspath, join
import sys

import numpy as np

from GoTools import Point, WriteG2
from GoTools.SurfaceFactory import LoftCurves
import GeoUtils.Interpolate as ip

from utils import grading, grading_double, gradspace
from airfoils import AirFoil


def fill_runner(af):
    af.fill()
    return af


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
        self.airfoils = [AirFoil.from_wd(self.p, wd) for wd in self.p.wingdef]
        self.p.dump_g2files('airfoils_raw', self.airfoils)


    def resample_airfoils(self):
        for af in self.airfoils:
            af.resample(self.p.n_te, self.p.n_back, self.p.n_front)
        self.p.dump_g2files('airfoils_resampled_radial', self.airfoils)


    def resolve_join(self):
        n = self.p.join_adds
        idx = self.p.join_index

        if n == 0:
            return

        for _ in xrange(n):
            new = AirFoil.average(self.airfoils[idx], self.airfoils[idx+1])
            self.airfoils.insert(idx + 1, new)

        for i in xrange(n):
            new = AirFoil.average(self.airfoils[idx+i-1], self.airfoils[idx+i])
            self.airfoils.insert(idx + i, new)

        del self.airfoils[idx + n]
        self.p.dump_g2files('airfoils_joined', self.airfoils)


    def resample_length(self):
        if self.p.length_mode == 'extruded' or self.p.mesh_mode == '2d':
            return

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
            self.airfoils.append(AirFoil.from_pts(pts, theta(z)))

        self.p.dump_g2files('airfoils_resampled_length', self.airfoils)


    def fill_airfoils(self):
        for af in self.airfoils:
            af.prepare_fill(self.p)

        progress('Filling slices', 0, len(self.airfoils))
        pool = mp.Pool(self.p.nprocs_mg)
        result = []
        for i, af in enumerate(pool.imap(fill_runner, self.airfoils)):
            result.append(af)
            progress('Filling slices', i+1, len(self.airfoils))

        self.airfoils = result
        self.p.dump_g2files('airfoils_filled', self.airfoils)


    def fill_sides(self):
        for af in self.airfoils:
            af.fill_sides()

        self.p.dump_g2files('airfoils_sides', self.airfoils)


    def extrude(self):
        assert(self.p.length_mode == 'extruded')

        zvals = np.linspace(0, self.p.length, self.p.n_length + 1)
        airfoil = self.airfoils[0]
        self.airfoils = [airfoil.translate(Point(0,0,z)) for z in zvals]


    def subdivide_airfoils(self):
        for af in self.airfoils:
            af.subdivide()

        self.p.dump_g2files('airfoils_subdivided', self.airfoils)


    def lower_order(self):
        for af in self.airfoils:
            af.lower_order(self.p.order)


    def output_semi3d(self):
        for i, af in enumerate(self.airfoils):
            path = abspath(join(self.p.out, 'slice-%03i' % (i+1)))
            af.output(path)

        zvals = [af.z() for af in self.airfoils]
        beam = ip.LinearCurve(pts=[Point(z,0,0) for z in zvals])
        WriteG2(join(self.p.out, 'beam.g2'), beam)

        self.p.out_yaml(join(self.p.out, 'parameters.yaml'))


    def output_2d(self):
        path = abspath(join(self.p.out, self.p.out))
        self.airfoils[0].output(path)


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
