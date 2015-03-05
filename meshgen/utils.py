import numpy as np

from itertools import chain
from math import pi, sin, cos
from numpy import matrix
from operator import *
import os

from GoTools import Point, Curve, Surface, Volume, WriteG2
from GoTools.VolumeFactory import LoftSurfaces, ExtrudeSurface
from GeoUtils.CurveUtils import GetCurvePoints
from GeoUtils.Elementary import Translate
from GeoUtils.Factory import LoftBetween
from GeoUtils.IO import InputFile
from GeoUtils.Refinement import UniformSurface
import GeoUtils.Interpolate as ip


ex = Point(1, 0, 0)
ey = Point(0, 1, 0)
ez = Point(0, 0, 1)


def convert_openfoam(path, delete=True):
    """Converts the IFEM xinp-file located at path (with no trailing extension) to OpenFOAM
    format. If delete is True, the IFEM files are removed.
    """
    f = InputFile('%s.xinp' % path)
    f.writeOpenFOAM(os.path.dirname(path))

    if delete:
        for postfix in ['.xinp', '.g2', '_nodenumbers.hdf5', '_nodenumbers.xml']:
            filename = path + postfix
            if os.path.exists(filename) and os.path.isfile(filename):
                os.remove(filename)


def extend_knots(knots):
    """Returns knots with the mean of the first and last two elements added."""
    return [knots[0], (knots[0]+knots[1])/2.] + knots[1:-1] + [(knots[-2]+knots[-1])/2., knots[-1]]


def standardize(obj):
    """Returns a new Surface of Curve with the same geometry but a standardized parametrization.
    The knots will be integers, i.e. [0, 1, 2, ...].
    """
    old_knots = obj.GetKnots() if type(obj) is Surface else [obj.GetKnots()]
    evl_knots = map(extend_knots, old_knots)
    param_knots = [extend_knots(range(len(kts))) for kts in old_knots]
    param_knots += [[0]*3 + range(len(kts)) + [len(kts)-1]*3 for kts in old_knots]

    pts = obj.EvaluateGrid(*evl_knots)
    if type(obj) is Curve:
        pts = matrix(pts)

    interpolator = ip.InterpolateCurve if type(obj) is Curve else ip.InterpolateSurface
    return interpolator(pts, *param_knots)


def merge_surfaces(srfa, srfb, dirc):
    """Merges two cubic surfaces with second order continuity along the direction given by dirc."""
    # The code will assume dirc == 0
    if dirc == 1:
        srfa.SwapParametrization()
        srfb.SwapParametrization()

    kua, kva = srfa.GetKnots(with_multiplicities=True)
    kub, _ = srfb.GetKnots(with_multiplicities=True)
    coeffsa = list(srfa)
    coeffsb = list(srfb)

    # Number of coefficients in the merging direction
    num = len(kva) - 4

    # Number of coefficients from each surface in the other direction
    numa, numb = len(kua) - 4, len(kub) - 4

    coeffs = []
    final = None
    for i in xrange(num):
        # Construct curves from each surface
        crva = Curve(4, kua, coeffsa[numa*i:numa*(i+1)], False)
        crvb = Curve(4, kub, coeffsb[numb*i:numb*(i+1)], False)

        # Merge them together to compute the coefficients
        crva.AppendCurve(crvb, continuity=2, reparam=False)
        crv = standardize(crva)
        coeffs += list(crv)

        # Grab the final knots only once
        if final is None:
            final = crv.GetKnots(with_multiplicities=True)

    # Form the new surface from the computed coefficients
    srf = Surface(4, 4, final, kva, coeffs, False)

    # Reset the direction
    if dirc == 1:
        srf.SwapParametrization()

    return srf


def mkcircle(center, radius, angle, nelems):
    """Generates a cubic nonrational approximation of a circle at the given center with the given
    radius. The parametrization starts at angle (given in degrees). The resulting curve will have
    nelems elements and is parametrized by curvelength.
    """
    alpha = angle * pi / 180
    thetas = np.linspace(0, 2*pi, nelems+1)
    pts = [center + Point(radius * cos(t + alpha),
                          radius * sin(t + alpha), 0)
           for t in thetas]

    return ip.CubicCurve(pts=pts, boundary=ip.PERIODIC).ReParametrize(0, 2*pi*radius)


def grading(length, ds, n, tol=1e-12, maxiters=1000):
    """Determines the grading factor required to fit n elements inside an interval of size length,
    using an initial element size of ds."""
    fr = lambda s, ds, n, r: s - ds*(1.0 - r**n) / (1.-r)
    dfrdr = lambda ds, n, r: ds * (r**(n-1) * (n*(1.-r) + r) - 1.0) / (1.-r)**2

    # Cell size for uniform mesh
    ds_unif = length / n

    # Find start guess for refinement factor
    r = 0.9 if ds_unif <= ds else 1.1

    # Newton loop
    its = 0
    eps = 10 * tol
    for it in xrange(1, maxiters + 1):
        if eps <= tol:
            break

        fval = fr(length, ds, n, r)
        r -= fval / dfrdr(ds, n, r)
        eps  = abs(fval)

    return r if it < maxiters else 0.0


def grading_double(length, d1, d2, N1, N2, tol=1e-12, maxiters=200):
    """Determines the grading factors required to fit N1 and N2 elements inside an interval of size
    length, where the middle elements are the same size, using initial element sizes d1 and d2 on
    both sides. The grading factors are towards the middle."""
    Ft = lambda r, d, N: d * (1.0 - r**N) / (1.0 - r)
    F = lambda r1, r2: Ft(r1, d1, N1) + Ft(r2, d2, N2) - length
    dFt = lambda r, d, N: (Ft(r, d, N) - d*N*r**(N-1)) / (1.0 - r)
    def dF(r1, r2):
        ret = dFt(r2, d2, N2)
        ret *= pow(r1**(N1-N2) * d1/d2, 1/(N2-1.0)) * (N1-1) / (N2-2)
        ret += dFt(r1, d1, N1)
        return ret

    ds_unif = length / (N1 + N2)
    r1 = 1.1 if ds_unif > d1 else 0.9
    r2 = lambda r1: pow(r1**(N1-1) * d1 / d2, 1/(N2-1.0))

    its = 0
    for _ in xrange(maxiters):
        f = F(r1, r2(r1))
        df = dF(r1, r2(r1))
        r1 -= f/df

        if abs(f) <= tol:
            break

        its += 1

    if its == maxiters:
        raise Exception("Maximal number of iterations reached")

    return r1, r2(r1)


def gradspace(start, step, factor, N):
    """Construct a graded point space from start, with the given initial step and grading factor.
    Returns N points."""
    def gen(start, step, factor, N):
        for _ in xrange(N):
            yield start
            start += step
            step *= factor
    return list(gen(float(start), float(step), float(factor), N))


def extend(edge, direction, distance, elements):
    """Generates a cubic surface between edge and edge translated along direction * distance.
    The resulting surface has the given number of elements in the new direction. Useful for
    extruding rectangles."""
    other = Translate(edge, direction * distance)
    surface = LoftBetween(edge, other)
    surface.RaiseOrder(0, 2)
    UniformSurface(surface, 2, elements-1)
    return surface


def is_patch(p):
    """Checks whether an object is a patch."""
    return type(p) in [Curve, Surface, Volume]


def is_list_of_patch(p):
    """Checks whether an object is a list of patches."""
    return type(p) is list and is_patch(p[0])


def deep(predicate, function):
    """Transforms the given function (which operates on some type determined by the predicate).
    The returned function instead operates on a nested list, preserving the list structure."""
    def deep_function(objs):
        if predicate(objs):
            return function(objs)
        else:
            return [deep_function(obj) for obj in objs]
    return deep_function


def deep_noret(predicate, function):
    """Performs the same operation as `deep`, except it discards returned values."""
    def deep_function(objs):
        if predicate(objs):
            function(objs)
        else:
            for obj in objs:
                deep_function(obj)
    return deep_function


def flatten_objects(patches):
    """Flattens a nested list of patches."""
    ret = []
    for p in patches:
        if is_patch(p):
            ret.append(p)
        elif type(p) is list:
            ret.extend(flatten_objects(p))
    return ret


def subdivide(patch, n, direction=0):
    """Subdivides a patch into n new patches along the given direction. Works for curves, surfaces
    and volumes. The new patches will be as close to uniformly sized as possible."""
    typ = type(patch)
    if typ is Curve:
        sub = lambda p, f, t: p.GetSubCurve(f[0], t[0])
    elif typ is Surface:
        sub = lambda p, f, t: p.GetSubSurf(f, t)
    elif typ is Volume:
        sub = lambda p, f, t: p.GetSubVol(f, t)

    if typ is Curve:
        nkts = len(patch.GetKnots())
    else:
        nkts = len(patch.GetKnots()[direction])
    indices = [i * (nkts - 1) / n for i in xrange(0, n+1)]

    kts = [patch.GetKnots()] if typ is Curve else patch.GetKnots()
    from_par = [k[0] for k in kts]
    to_par = [k[-1] for k in kts]
    parts = []
    for a, b in zip(indices[:-1], indices[1:]):
        from_par[direction] = kts[direction][a]
        to_par[direction] = kts[direction][b]
        parts.append(sub(patch, from_par, to_par))

    return parts


def deep_subdivide(patches, n, direction):
    """Nested list version of `subdivide`."""
    return deep(is_patch, lambda p: subdivide(p, n, direction))(patches)


def orient(patch, *args):
    """Orients a patch, performing the operations given by the additional arguments, in order. This
    function operates in-place. Valid operations are:
    - flip (for curves)
    - flipu, flipv (for surfaces and volumes)
    - flipw (for volumes)
    - swap (for surfaces)
    - swap[uvw][uvw] (for volumes)
    """
    for a in args:
        if a[:4] == 'flip':
            if type(patch) is Curve:
                patch.FlipParametrization()
            else:
                direction = 'uvw'.index(a[4])
                patch.FlipParametrization(direction)
        elif a == 'swap':
            if type(patch) is Surface:
                patch.SwapParametrization()
            else:
                dir1 = 'uvw'.index(a[4])
                dir2 = 'uvw'.index(a[5])
                patch.SwapParametrization(dir1, dir2)


def deep_orient(patches, *args):
    """Nested list version of `orient`."""
    deep_noret(is_patch, lambda p: orient(p, *args))(patches)


def lower_order(patch, target):
    """Lowers the orders of the given patch to the target order in all directions. Operates
    in-place."""
    orders = [patch.GetOrder()] if type(patch) is Curve else patch.GetOrder()
    lower = [order - target for order in orders]
    patch.LowerOrder(*lower)


def deep_lower_order(patches, target):
    """Nested list version of `lower_order`."""
    deep_noret(is_patch, lambda p: lower_order(p, target))(patches)


def extrude(patch):
    """Extrudes a surface with unit distance in the z-direction."""
    return ExtrudeSurface(patch, ez, 1.0)


def deep_extrude(patches):
    """Nested list version of `extrude`."""
    return deep(is_patch, extrude)(patches)


def deep_translate(patches, pt):
    """Nested list translation by a given point. Does not operate in-place."""
    return deep(is_patch, lambda p: Translate(p, pt))(patches)


def deep_index(patches, idx):
    """Indexes the last level(s) of a nested list structure of patches."""
    return deep(is_list_of_patch, itemgetter(idx))(patches)


def deep_loft(patches):
    """Lofts cubically the last level(s) of a nested list structure of surfaces."""
    if type(patches[0]) is Surface:
        return LoftSurfaces(patches, range(len(patches)), 4)
    else:
        lists = [list(q) for q in zip(*patches)]
        return [deep_loft(p) for p in lists]
