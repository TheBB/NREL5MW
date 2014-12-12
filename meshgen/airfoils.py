import numpy as np

from math import pi, sin, cos

from GoTools import Point, WriteG2
from GeoUtils.CurveUtils import CurveLengthParametrization, GetCurvePoints
import GeoUtils.Interpolate as ip
from GeoUtils.Refinement import UniformCurve

from utils import ez, mkcircle, grading_twosided, gradspace


def load_airfoil(filename, gap):
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

    def __init__(self, params, wd):
        if wd.foil == 'cylinder':
            center = Point(wd.chord * (.25 - wd.ao + wd.ac), 0, wd.z)
            self.curve = mkcircle(center, .5 * wd.chord, wd.theta, 200)
            self.len_te = params.len_te_cyl
            self.len_total = self.curve.GetKnots()[-1]
            self.len_upper = 0.5 * self.len_total
        else:
            self.init_normal(params, wd)


    def objects(self):
        return self.curve


    def init_normal(self, params, wd):
        xs, ys = load_airfoil(wd.foil, params.len_te / wd.chord)
        xs = wd.chord * (xs - (.25 + wd.ao - wd.ac))
        ys *= wd.chord

        pts_gap = [Point(x, y, wd.z).Rotate(ez, wd.theta * pi / 180) for x, y in zip(xs, ys)]
        tangents = [(pts_gap[1] - pts_gap[0]).Normalize(),
                    (pts_gap[-1] - pts_gap[-2]).Normalize()]
        pts_te = te_curve(pts_gap[0], pts_gap[-1], tangents[0], tangents[1])
        
        imid = (len(pts_te) - 1) / 2
        pts = pts_te[imid:-1] + pts_gap + pts_te[1:imid+1]
        self.curve = ip.CubicCurve(pts=pts, boundary=ip.PERIODIC)

        knots = self.curve.GetKnots()
        self.len_te = CurveLengthParametrization(pts_te)[-1]
        self.len_total = knots[-1]
        self.len_upper = knots[(len(knots) - 1) / 2]


    def resample(self, params):
        ds_back  = self.len_te / params.n_te / 2
        ds_front = ds_back

        r2, r3 = grading_twosided(self.len_upper - self.len_te / 2,
                                   ds_back, ds_front, params.n_back-1, params.n_front-1)
        r4, r5 = grading_twosided(self.len_total - self.len_upper - self.len_te / 2,
                                   ds_front, ds_back, params.n_front-1, params.n_back-1)

        knots = []
        last_ds = lambda: knots[-1] - knots[-2]
        last = lambda: knots[-1]

        knots += list(np.linspace(0, self.len_te / 2, params.n_te + 1))
        knots += gradspace(last(), ds_back, 1./r2, params.n_back  + 1)[1:]
        knots += gradspace(last(), last_ds(), r3,    params.n_front + 1)[1:]
        knots += gradspace(last(), last_ds(), 1./r4, params.n_front + 1)[1:]
        knots += gradspace(last(), last_ds(), r5,    params.n_back  + 1)[1:]
        knots += list(np.linspace(last(), self.len_total, params.n_te + 1))[1:]

        pts = [self.curve.Evaluate(k) for k in knots]
        self.curve = ip.CubicCurve(pts=pts, t=list(np.linspace(0, 1, len(pts))), boundary=ip.PERIODIC)

        length = CurveLengthParametrization(pts)
        self.len_te = length[params.n_te]
        self.len_upper = length[params.n_te + params.n_back + params.n_front]
        self.len_total = length[-1]
