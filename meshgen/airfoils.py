import numpy as np

from math import pi, sin, cos

from GoTools import Point, WriteG2
from GeoUtils.CurveUtils import CurveLengthParametrization, GetCurvePoints
import GeoUtils.Interpolate as ip
from GeoUtils.Refinement import UniformCurve


ez = Point(0, 0, 1)


def mkcircle(center, radius, angle, nelems):
    alpha = angle * pi / 180
    thetas = np.linspace(0, 2*pi, nelems+1)
    pts = [center + Point(radius * cos(t + alpha),
                          radius * sin(t + alpha), 0)
           for t in thetas]

    return ip.CubicCurve(pts=pts, boundary=ip.PERIODIC).ReParametrize(0, 2*pi*radius)


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

        WriteG2('out/test.g2', [self.curve])

    def objects(self):
        return self.curve
