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


def convert_openfoam(path, subfolder=False):
    f = InputFile('%s.xinp' % path)
    f.writeOpenFOAM(os.path.dirname(path))

    for postfix in ['.xinp', '.g2', '_nodenumbers.hdf5', '_nodenumbers.xml']:
        filename = path + postfix
        if os.path.exists(filename) and os.path.isfile(filename):
            os.remove(filename)


def extend_knots(knots):
    return [knots[0], (knots[0]+knots[1])/2.] + knots[1:-1] + [(knots[-2]+knots[-1])/2., knots[-1]]


def standardize(obj):
    old_knots = obj.GetKnots() if type(obj) is Surface else [obj.GetKnots()]
    evl_knots = map(extend_knots, old_knots)
    param_knots = [extend_knots(range(len(kts))) for kts in old_knots]
    param_knots += [[0]*3 + range(len(kts)) + [len(kts)-1]*3 for kts in old_knots]

    pts = obj.EvaluateGrid(*evl_knots)
    if type(obj) is Curve:
        pts = matrix(pts)

    interpolator = ip.InterpolateCurve if type(obj) is Curve else ip.InterpolateSurface
    return interpolator(pts, *param_knots)


def add_if_has(obj, attrib, lst):
    if hasattr(obj, attrib):
        lst.append(getattr(obj, attrib))


def merge_surfaces(srfa, srfb, dirc):
    if dirc == 1:
        srfa.SwapParametrization()
        srfb.SwapParametrization()

    kua, kva = srfa.GetKnots(with_multiplicities=True)
    kub, _ = srfb.GetKnots(with_multiplicities=True)
    coeffsa = list(srfa)
    coeffsb = list(srfb)

    num = len(kva) - 4
    numa, numb = len(kua) - 4, len(kub) - 4

    coeffs = []
    final = None
    for i in xrange(num):
        crva = Curve(4, kua, coeffsa[numa*i:numa*(i+1)], False)
        crvb = Curve(4, kub, coeffsb[numb*i:numb*(i+1)], False)

        crva.AppendCurve(crvb, continuity=2, reparam=False)
        crv = standardize(crva)
        coeffs += list(crv)

        if final is None:
            final = crv.GetKnots(with_multiplicities=True)

    srf = Surface(4, 4, final, kva, coeffs, False)

    if dirc == 1:
        srf.SwapParametrization()

    return srf


def mkcircle(center, radius, angle, nelems):
    alpha = angle * pi / 180
    thetas = np.linspace(0, 2*pi, nelems+1)
    pts = [center + Point(radius * cos(t + alpha),
                          radius * sin(t + alpha), 0)
           for t in thetas]

    return ip.CubicCurve(pts=pts, boundary=ip.PERIODIC).ReParametrize(0, 2*pi*radius)


def grading(length, ds, n, tol=1e-12, maxiters=1000):
    fr = lambda s, ds, n, r: s - ds*(1.0 - r**n) / (1.-r)
    dfrdr = lambda ds, n, r: ds * (r**(n-1) * (n*(1.-r) + r) - 1.0) / (1.-r)**2

    # cell size for uniform mesh
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


def _F(s, dx1, dx2, n1, n2, r10, r20):
    r2 = r20
    n1i = 1.0 / n1
    coeff1 = dx1 / dx2
    r2n2 = pow(r2, n2)
    r1 = pow(coeff1*r2n2, n1i)
    dx = dx2 / r2n2

    return s - dx * ((1.0-pow(r1,n1+1)) / (1.0-r1) + (1.0-pow(r2,n2+1)) / (1.0-r2))


def _dFdr2(dx1, dx2, n1, n2, r10, r20):
    r2 = r20
    r2n2 = pow(r2, n2)
    r2n2p1 = pow(r2, n2+1)
    c2 = (1.0-r2n2p1) / (1.0-r2)
    dc2dr2 = (-(n2+1) * r2n2 * (1.0-r2) - (1.0-r2n2p1)) / pow(1.0-r2,2)

    n1i = 1.0 / n1
    r1 = pow(dx1 / dx2 * r2n2, n1i)

    r1n1   = pow(r1, n1)
    r1n1p1 = pow(r1, n1+1)
    c1 = (1.0-r1n1p1) / (1.0-r1)
    dc1dr1 = (-(n1+1) * r1n1 * (1.0-r1) - (1.0-r1n1p1)) / pow(1.0-r1,2)
    dr1dr2 = n2 / n1 * pow(dx1/dx2, n1i) * pow(r2, n2/n1-1.0)

    dx = dx2 / r2n2
    dxdr2 = -n2 * dx2 / r2n2p1

    return - dxdr2*(c1 + c2) - dx*(dc1dr1*dr1dr2 + dc2dr2)


def grading_double(length, dx1, dx2, n1, n2, r1=1.1, r2=.9, tol=1e-12, maxiters=200):
    rf = lambda r2: pow(dx1 / dx2 * pow(r2, n2), 1./n1)

    eps = 10.0 * tol
    r1 = rf(r2)

    for it in xrange(1, maxiters + 1):
        if eps <= tol:
            break

        P = _F(length, dx1, dx2, n1, n2, r1, r2)
        dP = _dFdr2(dx1, dx2, n1, n2, r1, r2)
        r2 -= P / dP
        r1 = rf(r2)
        eps = abs(P)

    return (r1, r2) if it < maxiters else (0., 0.)


def gradspace(start, step, factor, N):
    def gen(start, step, factor, N):
        for _ in xrange(N):
            yield start
            start += step
            step *= factor
    return list(gen(float(start), float(step), float(factor), N))


def extend(edge, direction, distance, elements):
    other = Translate(edge, direction * distance)
    surface = LoftBetween(edge, other)
    surface.RaiseOrder(0, 2)
    UniformSurface(surface, 2, elements-1)
    return surface


is_patch = lambda p: type(p) in [Curve, Surface, Volume]
is_list_of_patch = lambda p: type(p) is list and is_patch(p[0])


def flatten_objects(patches):
    ret = []
    for p in patches:
        if is_patch(p):
            ret.append(p)
        elif type(p) is list:
            ret.extend(flatten_objects(p))
    return ret


def subdivide(patch, n, direction=0):
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


def orient(patch, *args):
    for a in args:
        if a[:4] == 'flip':
            direction = 'uvw'.index(a[4])
            patch.FlipParametrization(direction)
        elif a == 'swap':
            patch.SwapParametrization()


def lower_order(patch, target):
    orders = [order - target for order in patch.GetOrder()]
    patch.LowerOrder(*orders)


def extrude(patch):
    return ExtrudeSurface(patch, ez, 1.0)


def deep(predicate, function):
    def deep_function(objs):
        if predicate(objs):
            return function(objs)
        else:
            return [deep_function(obj) for obj in objs]
    return deep_function


def deep_noret(predicate, function):
    def deep_function(objs):
        if predicate(objs):
            function(objs)
        else:
            for obj in objs:
                deep_function(obj)
    return deep_function


def deep_subdivide(patches, n, direction):
    return deep(is_patch, lambda p: subdivide(p, n, direction))(patches)


def deep_orient(patches, *args):
    deep_noret(is_patch, lambda p: orient(p, *args))(patches)


def deep_translate(patches, pt):
    return deep(is_patch, lambda p: Translate(p, pt))(patches)


def deep_index(patches, idx):
    return deep(is_list_of_patch, itemgetter(idx))(patches)


def deep_lower_order(patches, target):
    deep_noret(is_patch, lambda p: lower_order(p, target))(patches)


def deep_extrude(patches):
    return deep(is_patch, extrude)(patches)


def deep_loft(patches):
    if type(patches[0]) is Surface:
        return LoftSurfaces(patches, range(len(patches)), 4)
    else:
        lists = [list(q) for q in zip(*patches)]
        return [deep_loft(p) for p in lists]
