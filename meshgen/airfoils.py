import numpy as np

from math import pi, sin, cos

from GoTools import Point, WriteG2
from GeoUtils.CurveUtils import CurveLengthParametrization, GetCurvePoints
import GeoUtils.Interpolate as ip
from GeoUtils.Refinement import UniformCurve

from utils import ez, mkcircle, grading_double, gradspace


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

    def __init__(self):
        pass


    @classmethod
    def from_wd(params, wd):
        obj = cls()
        obj.theta = wd.theta

        if wd.foil == 'cylinder':
            center = Point(wd.chord * (.25 - wd.ao + wd.ac), 0, wd.z)
            obj.curve = mkcircle(center, .5 * wd.chord, wd.theta, 200)
            obj.len_te = params.len_te_cyl
            obj.len_total = self.curve.GetKnots()[-1]
            obj.len_upper = 0.5 * self.len_total
        else:
            obj.init_normal(params, wd)
        


    @classmethod
    def from_pts(cls, pts, theta):
        obj = cls()
        obj._from_pts(pts)
        obj.theta = theta


    @classmethod
    def average(cls, afa, afb):
        assert(afa.curve.GetKnots() == afb.curve.GetKnots())

        pts = [(pta + ptb) / 2 for pta, ptb in zip(GetCurvePoints(afa.curve),
                                                   GetCurvePoints(afb.curve))]
        return cls(pts=pts, theta=(afa.theta + afb.theta) / 2)


    def objects(self):
        return self.curve


    def z(self):
        return self.curve.Evaluate(self.curve.GetKnots()[0])[2]


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

        self.len_te = CurveLengthParametrization(pts_te)[-1]


    def resample(self, n_te, n_back, n_front):
        knots = self.curve.GetKnots()
        len_total = knots[-1]
        len_upper = knots[(len(knots) - 1) / 2]

        ds_back  = self.len_te / n_te / 2
        ds_front = ds_back

        r2, r3 = grading_double(len_upper - self.len_te / 2,
                                ds_back, ds_front, n_back-1, n_front-1)
        r4, r5 = grading_double(len_total - len_upper - self.len_te / 2,
                                ds_front, ds_back, n_front-1, n_back-1)

        knots = []
        last_ds = lambda: knots[-1] - knots[-2]
        last = lambda: knots[-1]

        knots += list(np.linspace(0, self.len_te / 2, n_te + 1))
        knots += gradspace(last(), last_ds(), 1./r2, n_back  + 1)[1:]
        knots += gradspace(last(), last_ds(), r3,    n_front + 1)[1:]
        knots += gradspace(last(), last_ds(), 1./r4, n_front + 1)[1:]
        knots += gradspace(last(), last_ds(), r5,    n_back  + 1)[1:]
        knots += list(np.linspace(last(), len_total, n_te + 1))[1:]

        pts = [self.curve.Evaluate(k) for k in knots]
        self._from_pts(pts)

        del self.len_te


    def prepare_tfi(self, radius):
        x, y, _ = self.curve.Evaluate(0)
        print x, y

        theta = np.arctan(y / x)
        print 180 * theta / pi
        # self.circle = mkcircle(Point(0, 0, self.z()),
        #                        radius, self.theta/2,
        #                        len(self.curve.GetKnots()) - 1)
        # self.curve = mkcircle(center, .5 * wd.chord, wd.theta, 200)


    # def tfi(self):
    #     pts = [Point(x,0,0) for x in np.linspace(0,1,10000)]
    #     return ip.LinearCurve(pts=pts)


    def _from_pts(self, pts):
        params = list(np.linspace(0, 1, len(pts)))
        self.curve = ip.CubicCurve(pts=pts, t=params, boundary=ip.PERIODIC)
