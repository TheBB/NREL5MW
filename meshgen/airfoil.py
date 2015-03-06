import os

import numpy as np

from GoTools import Point
from GeoUtils.CurveUtils import CurveLengthParametrization, GetCurvePoints
from GeoUtils.Refinement import UniformCurve
import GeoUtils.Interpolate as ip

from utils import *


def load_airfoil(filename, gap):
    """Loads an airfoil from the given file (relative to the airfoils folder).  Introduces a
    trailing edge gap of the given size.  Returns x-values and y-values."""
    filename = os.path.join(os.path.dirname(__file__), '../airfoils', filename)
    data = np.loadtxt(filename)

    # Determine the x- and y-values of the upper and lower sides
    angles = np.arctan(data[:,3])
    x_up = data[:,0] - .5 * data[:,1] * np.sin(angles)
    y_up = data[:,2] + .5 * data[:,1] * np.cos(angles)
    x_dn = data[:,0] + .5 * data[:,1] * np.sin(angles)
    y_dn = data[:,2] - .5 * data[:,1] * np.cos(angles)

    # Determine the midpoints and a normalized curve length parametrization
    midpoints = [Point(x, y, 0) for x, y in zip((x_up + x_dn)/2, (y_up + y_dn)/2)]
    coeffs = np.array(CurveLengthParametrization(midpoints, normalize=True))

    # Modify to produce the right gap
    diff = (y_up[-1] - y_dn[-1] - gap) / 2
    y_up -= coeffs * diff
    y_dn += coeffs * diff

    # Return points in positive direction from the trailing edge
    xs = np.hstack((np.flipud(x_up), x_dn[1:]))
    ys = np.hstack((np.flipud(y_up), y_dn[1:]))
    return xs, ys


def te_curve(pta, ptb, tnga, tngb):
    """Generates a rounded trailing edge between pta (above) and ptb (below), given the tangents of
    the airfoil at those points (tnga pointing 'forward' and tngb pointing 'backward'). Returns a
    list of points describing the rounded trailing edge, starting with ptb and ending with pta."""
    gap = abs(pta - ptb)

    # More or less the backward pointing normal
    norm = (tngb - tnga).Normalize()

    # The middle point of the curve
    ptm = (pta + ptb + gap * norm) / 2

    curve = ip.CubicCurve(pts=[ptb, ptm, pta], boundary=ip.TANGENT, der=[tngb, tnga])
    UniformCurve(curve, n=5)

    return GetCurvePoints(curve)


class AirFoil(object):
    """This class represents a single one-dimensional periodic airfoil curve."""

    def __init__(self, params, theta):
        """Private constructor. Foreign code should call one of the specialized constructors
        (`from_wingdef`, `from_pts` or `from_mean`)."""
        self.p = params
        self.theta = theta

    @classmethod
    def from_wingdef(cls, params, wd):
        """Produces an airfoil from a wing definition."""
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
        """Produces an airfoil from a point cloud (starting and ending at the trailing edge,
        positive rotation direction."""
        airfoil = cls(params, theta)
        airfoil.curve = ip.CubicCurve(pts=pts, t=range(len(pts)), boundary=ip.PERIODIC)
        return airfoil

    @classmethod
    def from_mean(cls, a, b):
        """Produces an airfoil by taking the mean of two other airfoils."""
        pts_a, pts_b = map(GetCurvePoints, [a.curve, b.curve])
        pts = [(p_a + p_b)/2 for p_a, p_b in zip(pts_a, pts_b)]
        return cls.from_pts(a.p, (a.theta + b.theta)/2, pts)

    def translate(self, pt):
        """Creates a new airfoil by translating this one by a given vector."""
        airfoil = AirFoil(self.p, self.theta)
        airfoil.curve = deep_translate(self.curve, pt)
        return airfoil

    def objects(self):
        """Required for debug output."""
        return [self.curve]

    def z(self):
        """Returns the z-location of this airfoil."""
        return self.curve[0][2]

    def resample(self):
        """Resamples the airfoil according to the given parameters. The airfoil must have been
        constructed using `from_wingdef`. It is an error to call `resample` twice."""
        # The number of elements in the trailing edge, front and back parts
        n_te, n_back, n_front = self.p.n_te, self.p.n_back, self.p.n_front

        # Get the lengths of the different parts
        knots = self.curve.GetKnots()
        len_total = knots[-1]
        len_upper = knots[(len(knots) - 1) / 2]

        # Element sizes 
        ds_back  = self.len_te / n_te / 2
        ds_front = len_total / (n_te + n_back + n_front) / 2 * 1e-1

        # Grading factors
        r2, r3 = grading_double(len_upper - self.len_te / 2,
                                ds_back, ds_front, n_back, n_front)
        r4, r5 = grading_double(len_total - len_upper - self.len_te / 2,
                                ds_front, ds_back, n_front, n_back)

        # Construct the knot vector used for evaluation
        knots = []
        last_ds = lambda: knots[-1] - knots[-2]
        last = lambda: knots[-1]

        knots += list(np.linspace(0, self.len_te / 2, n_te + 1))
        knots += gradspace(last(), last_ds(), r2,    n_back  + 1)[1:]
        knots += gradspace(last(), last_ds(), 1./r3, n_front + 1)[1:]
        knots += gradspace(last(), last_ds(), r4,    n_front + 1)[1:]
        knots += gradspace(last(), last_ds(), 1./r5, n_back  + 1)[1:]
        knots += list(np.linspace(last(), len_total, n_te + 1))[1:]

        # Evaluate and interpolate with new knot vector
        pts = [self.curve.Evaluate(k) for k in knots]
        self.curve = ip.CubicCurve(pts=pts, t=range(len(pts)), boundary=ip.PERIODIC)

        # `len_te` is no longer needed, delete it to produce an error if `resample` is called again
        del self.len_te
