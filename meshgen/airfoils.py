import numpy as np

from math import pi, sin, cos, sqrt

from GoTools import Point, WriteG2, Curve
from GoTools.CurveFactory import Circle, IntersectCurve, LineSegment, NonRationalCurve
from GeoUtils.CurveUtils import CurveLengthParametrization, GetCurvePoints
from GeoUtils.Factory import LoftBetween
from GeoUtils.Refinement import UniformCurve, GeometricRefineCurve, GeometricRefineSurface
import GeoUtils.Interpolate as ip
import GeoUtils.TFI as tfi

from utils import ex, ez, merge_surfaces, mkcircle, grading, grading_double, gradspace


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
    def from_wd(cls, params, wd):
        obj = cls()
        obj.theta = wd.theta

        if wd.foil == 'cylinder':
            center = Point(wd.chord * (.25 - wd.ao + wd.ac), 0, wd.z)
            obj.curve = mkcircle(center, .5 * wd.chord, wd.theta, 200)
            obj.len_te = params.len_te_cyl
            obj.len_total = obj.curve.GetKnots()[-1]
            obj.len_upper = 0.5 * obj.len_total
        else:
            obj.init_normal(params, wd)

        return obj


    @classmethod
    def from_pts(cls, pts, theta):
        obj = cls()
        obj._from_pts(pts)
        obj.theta = theta

        return obj


    @classmethod
    def average(cls, afa, afb):
        assert(afa.curve.GetKnots() == afb.curve.GetKnots())

        pts = [(pta + ptb) / 2 for pta, ptb in zip(GetCurvePoints(afa.curve),
                                                   GetCurvePoints(afb.curve))]
        return cls.from_pts(pts, (afa.theta + afb.theta) / 2)


    def objects(self):
        objects = []

        if hasattr(self, 'split'):
            objects += self.split
        else:
            objects.append(self.curve)

        return objects


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
        ds_front = len_total / (n_te + n_back + n_front) / 2 * 1e-1

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


    def prepare_fill(self, params):
        self.p = params


    def fill(self):
        trailing = self._make_trailing()
        inner = self._make_inner(trailing)
        split_inner = self._split_inner(inner)
        split_middle = self._make_middle(split_inner)
        split = self._merge(split_inner, split_middle)

        self.split = split


    def _make_trailing(self):
        k = self.curve.GetKnots()[0]
        p_inner = self.curve.Evaluate(k)
        n_inner = self.curve.EvaluateTangent(k).Rotate(ez, -pi/2).Normalize()

        p_outer = Point(2 * self.p.radius, 0, self.z())

        trailing = ip.CubicCurve(pts=[p_inner, p_outer], der=[n_inner, ex],
                                 boundary=ip.TANGENT)

        circle = Circle(ez * self.z(), self.p.radius, ez)
        point = IntersectCurve(trailing, circle)[1][0]
        knot = trailing.GetParameterAtPoint(point)[0]

        trailing.InsertKnot(knot)
        trailing = trailing.GetSubCurve(0, knot)
        fac = self.p.radial_grading(knot)
        GeometricRefineCurve(trailing, fac, self.p.n_circle - 1)

        self.inner_fac = fac

        return trailing


    def _make_inner(self, trailing):
        point = trailing.Evaluate(trailing.GetKnots()[-1])
        theta = np.arctan(point[1] / point[0]) * 180 / pi
        radius = abs(trailing.Evaluate(trailing.GetKnots()[-1]) - ez * self.z())
        
        circle =  mkcircle(ez * self.z(), radius, theta,
                           len(self.curve.GetKnots()) - 1)

        normalf = lambda _, crv, kts: [crv.EvaluateTangent(k)
                                          .Rotate(ez, -pi/2)
                                          .Normalize()
                                       for k in kts]
        init_normalf = lambda crv, kts: normalf(None, crv, kts)

        kwargs = {}
        if self.p.smoothing:
            kwargs = {'fac_smooth': .98,
                      'nsweeps': 1,
                      'ranges': [(0,2*self.p.n_te), (-2*self.p.n_te,0)]}

        return tfi.OrthogonalSurface(
            [self.curve, circle, trailing, trailing],
            init_normalf, normalf, fac_blend=.9, bnd_layer=self.p.n_bndlayer,
            **kwargs
        )


    def _split_inner(self, inner):
        ksu, ksv = inner.GetKnots()
        ksu_split = [ksu[i] for i in self.p.angular_splits()]

        splits = []
        for kp, kn in zip(ksu_split[:-1], ksu_split[1:]):
            splits.append(inner.GetSubSurf((kp, ksv[0]), (kn, ksv[-1])))

        return splits


    def _make_middle(self, inners):
        radius = 2 * self.p.radius
        z = self.z()

        outer_pts = [Point(radius, 0, z),
                     Point(radius, radius, z),
                     Point(0, radius, z),
                     Point(-radius, radius, z),
                     Point(-radius, 0, z),
                     Point(-radius, -radius, z),
                     Point(0, -radius, z),
                     Point(radius, -radius, z),
                     Point(radius, 0, z)]

        kus, kvs = inners[0].GetKnots()
        ipt = inners[0].Evaluate(kus[0], kvs[-1])
        length = abs(ipt - outer_pts[0])
        ds = abs(ipt - inners[0].Evaluate(kus[0], kvs[-2])) * self.inner_fac
        fac = grading(length, ds, self.p.n_square)

        outers = []
        for opp, opn, inner in zip(outer_pts[:-1], outer_pts[1:], inners):
            kus, kvs = inner.GetKnots()

            temp_curve = LineSegment(opp, opn)
            UniformCurve(temp_curve, len(kus) - 2)
            out = ip.CubicCurve(pts=GetCurvePoints(temp_curve), boundary=ip.NATURAL, t=kus)

            surface = LoftBetween(inner.GetEdges()[2], out).RaiseOrder(0, 2)
            GeometricRefineSurface(surface, 2, fac, self.p.n_square - 1)

            diff = surface.GetKnots()[1][1]
            surface.ReParametrize(kus[0], kus[-1], kvs[-1], kvs[-1] + 1/diff)

            outers.append(surface)

        return outers


    def _merge(self, inner, middle):
        return [merge_surfaces(i, m, 1) for i, m in zip(inner, middle)]


    def _from_pts(self, pts):
        params = list(np.linspace(0, 1, len(pts)))
        self.curve = ip.CubicCurve(pts=pts, t=params, boundary=ip.PERIODIC)
