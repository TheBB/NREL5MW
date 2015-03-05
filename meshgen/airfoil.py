import os

import numpy as np

from GoTools import Point
from GeoUtils.CurveUtils import CurveLengthParametrization, GetCurvePoints
from GeoUtils.Refinement import UniformCurve
import GeoUtils.Interpolate as ip

from utils import *


def load_airfoil(filename, gap):
    filename = os.path.join(os.path.dirname(__file__), '../airfoils', filename)
    data = np.loadtxt(filename)

    angles = np.arctan(data[:,3])
    x_up = data[:,0] - .5 * data[:,1] * np.sin(angles)
    y_up = data[:,2] + .5 * data[:,1] * np.cos(angles)
    x_dn = data[:,0] + .5 * data[:,1] * np.sin(angles)
    y_dn = data[:,2] - .5 * data[:,1] * np.cos(angles)

    midpoints = [Point(x, y, 0) for x, y in zip((x_up + x_dn)/2, (y_up + y_dn)/2)]
    coeffs = np.array(CurveLengthParametrization(midpoints, normalize=True))

    diff = (y_up[-1] - y_dn[-1] - gap) / 2
    y_up -= coeffs * diff
    y_dn += coeffs * diff

    xs = np.hstack((np.flipud(x_up), x_dn[1:]))
    ys = np.hstack((np.flipud(y_up), y_dn[1:]))
    return xs, ys


def te_curve(pta, ptb, tnga, tngb):
    gap = abs(pta - ptb)
    norm = (tngb - tnga).Normalize()
    ptm = (pta + ptb + gap * norm) / 2

    curve = ip.CubicCurve(pts=[ptb, ptm, pta], boundary=ip.TANGENT, der=[tngb, tnga])
    UniformCurve(curve, n=5)

    return GetCurvePoints(curve)


class AirFoil(object):

    def __init__(self, params, theta):
        self.p = params
        self.theta = theta

    @classmethod
    def from_wingdef(cls, params, wd):
        airfoil = cls(params, wd.theta)

        if wd.foil == 'cylinder':
            center = Point(wd.chord * (.25 - wd.ao + wd.ac), 0, wd.z)
            airfoil.curve = mkcircle(center, .5 * wd.chord, wd.theta, 200)
            airfoil.len_te = params.len_te_cyl
        else:
            xs, ys = load_airfoil(wd.foil, params.len_te / wd.chord)
            xs = wd.chord * (xs - (.25 + wd.ao - wd.ac))
            ys *= wd.chord

            pts_gap = [Point(x, y, wd.z).Rotate(ez, wd.theta * pi / 180) for x, y in zip(xs, ys)]
            tangents = [(pts_gap[1] - pts_gap[0]).Normalize(),
                        (pts_gap[-1] - pts_gap[-2]).Normalize()]
            pts_te = te_curve(pts_gap[0], pts_gap[-1], tangents[0], tangents[1])

            imid = (len(pts_te) - 1) / 2
            pts = pts_te[imid:-1] + pts_gap + pts_te[1:imid+1]

            airfoil.curve = ip.CubicCurve(pts=pts, boundary=ip.PERIODIC)
            airfoil.len_te = CurveLengthParametrization(pts_te)[-1]

        return airfoil

    @classmethod
    def from_pts(cls, params, theta, pts):
        airfoil = cls(params, theta)
        airfoil.curve = ip.CubicCurve(pts=pts, t=range(len(pts)), boundary=ip.PERIODIC)
        return airfoil

    @classmethod
    def from_mean(cls, a, b):
        pts_a, pts_b = map(GetCurvePoints, [a.curve, b.curve])
        pts = [(p_a + p_b)/2 for p_a, p_b in zip(pts_a, pts_b)]
        return cls.from_pts(a.p, (a.theta + b.theta)/2, pts)

    def translate(self, pt):
        airfoil = AirFoil(self.p, self.theta)
        airfoil.curve = deep_translate(self.curve, pt)
        return airfoil

    def objects(self):
        return [self.curve]

    def z(self):
        return self.curve[0][2]

    def resample(self):
        n_te, n_back, n_front = self.p.n_te, self.p.n_back, self.p.n_front

        knots = self.curve.GetKnots()
        len_total = knots[-1]
        len_upper = knots[(len(knots) - 1) / 2]

        ds_back  = self.len_te / n_te / 2
        ds_front = len_total / (n_te + n_back + n_front) / 2 * 1e-1

        r2, r3 = grading_double(len_upper - self.len_te / 2,
                                ds_back, ds_front, n_back, n_front)
        r4, r5 = grading_double(len_total - len_upper - self.len_te / 2,
                                ds_front, ds_back, n_front, n_back)

        knots = []
        last_ds = lambda: knots[-1] - knots[-2]
        last = lambda: knots[-1]

        knots += list(np.linspace(0, self.len_te / 2, n_te + 1))
        knots += gradspace(last(), last_ds(), r2,    n_back  + 1)[1:]
        knots += gradspace(last(), last_ds(), 1./r3, n_front + 1)[1:]
        knots += gradspace(last(), last_ds(), r4,    n_front + 1)[1:]
        knots += gradspace(last(), last_ds(), 1./r5, n_back  + 1)[1:]
        knots += list(np.linspace(last(), len_total, n_te + 1))[1:]

        pts = [self.curve.Evaluate(k) for k in knots]
        self.curve = ip.CubicCurve(pts=pts, t=range(len(pts)), boundary=ip.PERIODIC)

        del self.len_te
