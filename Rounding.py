from itertools import chain, cycle
from numpy import linspace, pi, sqrt, sin, arccos, angle, arctan, log, exp, array, matrix
from scipy import optimize

import sys

from GoTools import *
from GoTools.CurveFactory import *
from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from GeoUtils.Refinement import *

import GeoUtils.Interpolate as ip
import GeoUtils.TFI as tfi


def fr(s, ds, N, r):
    return s - ds*(1.0 - pow(r,N)) / (1.0-r)


def dfrdr(ds, N, r):
    value = pow(r,N-1) * (N*(1.0-r) + r) - 1.0
    return ds * value / pow(1.0-r,2)


def ComputeFactor(s, ds, N):
    # cell size for uniform mesh
    ds_unif = s / N

    # Find start guess for refinement factor
    r = 0.9 if ds_unif <= ds else 1.1

    # Newton loop
    maxIts = 1000
    tol    = 1.0e-10
    its    = 0
    eps    = 1.0

    while eps > tol and its < maxIts:
        fval = fr(s, ds, N, r)
        eps  = abs(fval)

        dfval = dfrdr(ds, N, r)
        if not dfval: return r
        dr    = -fval/dfval
        r     = r + dr

        its += 1

    return r if its < maxIts else 0.0


def CurveLengthParametrization(pts, normalize=False, xyplane=False):
    # Curve length parametrization
    s = 0.0
    knots = [s]
    for pp, pn in zip(pts[:-1], pts[1:]):
        if xyplane:
            s += abs(pn - pp - Point(0,0,pn[2]-pp[2]))
        else:
            s += abs(pn - pp)
        knots.append(s)

    # Normalized curve length parametrization
    if normalize:
        s = knots[-1]
        knots = [k/s for k in knots]

    return knots


def GetCurvePoints(curve):
    return [curve.Evaluate(k) for k in curve.GetKnots()]


def UniformParametrization(pts, normalize=False):
    ret = range(0, len(pts))
    if normalize:
        ret = [float(r)/ret[-1] for r in ret]
    return ret


def GradedSpace(start, step, factor, N):
    def GeneratorGradedSpace(start, step, factor, N):
        for _ in xrange(N):
            yield start
            start += step
            step *= factor
    return list(GeneratorGradedSpace(float(start), float(step), float(factor), N))


def CSeg(center, start, angle, axis, n, space='lin'):
    points = []

    if space == 'lin':
        space = linspace(0, angle, n+1)
    else:
        space = GradedSpace(0.0, angle * (1.0 - space) / (1.0 - space**n), space, n+1)

    space[-1] = angle

    for v in space:
        q = start.Clone()
        q = q - center
        q.Rotate(axis, v)
        q = q + center
        points.append(q)

    crv = InterpolateCurve(points, CurveLengthParametrization(points),
                           ((center - points[0]) % axis).Normalize(),
                           ((center - points[-1]) % axis).Normalize()) 

    return crv


def MkTopCurve(top, btm, topd, iparams):
    iparams = [(int(p*len(top)), h) for p,h in iparams]
    ipts = [top[0]] + [(top[i] + btm[i])/2 + Point(0,0,h) for i,h in iparams] + [top[-1]]

    topcurve = ip.CubicCurve(pts=ipts, t=list(linspace(0, 1, len(ipts))),
                             boundary=ip.TANGENT, der=[topd[0], -topd[-1]])

    ipts = [topcurve.Evaluate(k) for k in linspace(0, 1, 3*len(top))]
    knots = CurveLengthParametrization(ipts, True, True)
    topcurve = ip.CubicCurve(pts=ipts, t=knots, boundary=ip.TANGENT, der=[topd[0], -topd[-1]])

    ipts = [(t+b)/2 for t,b in zip(top, btm)]
    params = CurveLengthParametrization(ipts, True, True)
    ipts = [ipt + Point(0,0,topcurve.Evaluate(k)[2]-ipt[2]) for k, ipt in zip(params, ipts)]
    lens = CurveLengthParametrization(ipts, True)
    return ip.CubicCurve(pts=ipts, t=params, boundary=ip.TANGENT, der=[topd[0], -topd[-1]]), params, lens


def ProjectSurfaceToSurface(source, target, ulow, uhigh):

    tkvs = target.GetKnots()[1]

    def ProjectPointToSurface(p, surface):
        u, v = 0.0, 0.0
        udir = True

        while True:
            if udir:
                f = lambda x: abs(target.Evaluate(x, v) - p)
                u = optimize.fminbound(f, ulow, uhigh)
            else:
                f = lambda x: abs(target.Evaluate(u, x) - p)
                vt = optimize.fminbound(f, tkvs[0], tkvs[-1])
                if abs(v - vt) < 1e-4:
                    break
                v = vt

            udir = not udir

        return target.Evaluate(u, v)

    kus, kvs = source.GetKnots()
    ekus = [kus[0], (kus[0]+kus[1])/2] + kus[1:-1] + [(kus[-2]+kus[-1])/2, kus[-1]]
    ekvs = [kvs[0], (kvs[0]+kvs[1])/2] + kvs[1:-1] + [(kvs[-2]+kvs[-1])/2, kvs[-1]]
    pts = [ProjectPointToSurface(source.Evaluate(ku, kv), target) for kv in ekvs for ku in ekus]

    return ip.InterpolateSurface(pts, ekus, ekvs,
                                 [kus[0]]*3 + kus + [kus[-1]]*3,
                                 [kvs[0]]*3 + kvs + [kvs[-1]]*3).ReParametrize()


# I have made a drawing that turns this incomprehensible pile of garbage
# into a nirvana of enlightenment.  Consult me for salvation.  Salvation
# not guaranteed.  -- Eivind
def MkSurfaces(N1, N2, params, distfun, hmSurf, hmCurve, height):
    def ipl(c, ps):
        c.ReParametrize()
        eps = [ps[0], (ps[0]+ps[1])/2] + ps[1:-1] + [(ps[-2]+ps[-1])/2, ps[-1]]
        pts = [c.Evaluate(p) for p in eps]
        kts = [ps[0]]*3 + ps + [ps[-1]]*3
        return ip.InterpolateCurve(matrix(pts), eps, kts)

    def clparam(c, ps):
        c = ipl(c, ps)

        hcurve = hmCurve(c)
        hkts = hcurve.GetKnots()
        hkts = [hkts[0], (hkts[0]+hkts[1])/2] + hkts[1:-1] + [(hkts[-2]+hkts[-1])/2, hkts[-1]]
        nkts = CurveLengthParametrization([hcurve.Evaluate(k) for k in hkts], normalize=True)
        tkts = [nkts[0]]*4 + nkts[2:-2] + [nkts[-1]]*4

        curve = ip.InterpolateCurve([c.Evaluate(k) for k in hkts], nkts, tkts)
        return ipl(curve, ps)

    def flip(c):
        return c.Clone().FlipParametrization()

    N = 2*N1 + 2*N2
    p1 = list(linspace(0, 1, N1+1))
    p2 = list(linspace(0, 1, N2+1))

    a = 0.0
    b = float(N1)
    c = float(N1 + N2)
    d = float(N1 + 2*N2)
    e = float(2*N1 + 2*N2)


    temp = hmCurve(clparam(LineSegment(Point(1.0, c, 0), Point(0.0, c, 0)),
                           list(linspace(0, 1, 100))))
    dist = CurveLengthParametrization(GetCurvePoints(temp))[-1]
    N0 = int(ceil(dist / params.dzTip)) - 2*N2
    p0 = list(linspace(0, 1, N0+1))

    n0 = float(N0) / (N0+N2) / 2 * N
    n2 = float(N0 + 2*N2) / (N0+N2) / 2 * N
    hl = 0.5 * float(N0) / (N0 + N1)

    # slide and fac are magic numbers, feel free to experiment
    # They determine the location of the center point in each quadrant
    slide = float(N0) / (N0 + N2)
    fac = 1 - slide/3

    mb = fac * ((1.0 - slide) * b + slide * n0) + (1.0 - fac) * c
    md = fac * ((1.0 - slide) * d + slide * n2) + (1.0 - fac) * c

    th = fac * pi/2

    ml1 = clparam(LineSegment(Point(1.0, c, 0), Point(1.0 - hl, c, 0)), p0)
    ml2 = ipl(LineSegment(Point(1.0 - hl, c, 0), Point(0.5, c, 0)), p1)
    mr1 = clparam(LineSegment(Point(0.0, c, 0), Point(hl, c, 0)), p0)
    mr2 = ipl(LineSegment(Point(hl, c, 0), Point(0.5, c, 0)), p1)

    t0 = clparam(LineSegment(Point(0.5, a, 0), Point(0.5, n0, 0)), p0)
    t1 = ipl(LineSegment(Point(0.5, n0, 0), Point(0.5, c, 0)), p2)
    t2 = ipl(LineSegment(Point(0.5, n2, 0), Point(0.5, c, 0)), p2)
    t3 = clparam(LineSegment(Point(0.5, e, 0), Point(0.5, n2, 0)), p0)

    lt0 = clparam(ip.CubicCurve(pts=[Point(1.0, b, 0), Point(1-hl, mb, 0)],
                                t=[0, 1], boundary=ip.TANGENT,
                                der=[slide * Point(-1, 0, 0),
                                     0.01 * Point(-sin(th), cos(th)*N, 0)]), p0)
    lt1 = ipl(LineSegment(Point(1-hl, mb, 0), Point(1-hl, c, 0)), p2)
    lt2 = ipl(LineSegment(Point(1-hl, md, 0), Point(1-hl, c, 0)), p2)
    lt3 = clparam(ip.CubicCurve(pts=[Point(1.0, d, 0), Point(1-hl, md, 0)],
                                t=[0, 1], boundary=ip.TANGENT,
                                der=[slide * Point(-1, 0, 0),
                                     0.01 * Point(-sin(th), -cos(th)*N, 0)]), p0)

    rt0 = clparam(ip.CubicCurve(pts=[Point(0.0, b, 0), Point(hl, mb, 0)],
                                t=[0, 1], boundary=ip.TANGENT,
                                der=[slide * Point(1, 0, 0),
                                     0.01 * Point(sin(th), cos(th)*N, 0)]), p0)
    rt1 = ipl(LineSegment(Point(hl, mb, 0), Point(hl, c, 0)), p2)
    rt2 = ipl(LineSegment(Point(hl, md, 0), Point(hl, c, 0)), p2)
    rt3 = clparam(ip.CubicCurve(pts=[Point(0.0, d, 0), Point(hl, md, 0)],
                                t=[0, 1], boundary=ip.TANGENT,
                                der=[slide * Point(1, 0, 0),
                                     0.01 * Point(sin(th), -cos(th)*N, 0)]), p0)

    r0 = ipl(LineSegment(Point(0.0, a, 0), Point(0.0, b, 0)), p1)
    r1 = ipl(LineSegment(Point(0.0, b, 0), Point(0.0, c, 0)), p2)
    r2 = ipl(LineSegment(Point(0.0, d, 0), Point(0.0, c, 0)), p2)
    r3 = ipl(LineSegment(Point(0.0, e, 0), Point(0.0, d, 0)), p1)

    l0 = ipl(LineSegment(Point(1.0, a, 0), Point(1.0, b, 0)), p1)
    l1 = ipl(LineSegment(Point(1.0, b, 0), Point(1.0, c, 0)), p2)
    l2 = ipl(LineSegment(Point(1.0, d, 0), Point(1.0, c, 0)), p2)
    l3 = ipl(LineSegment(Point(1.0, e, 0), Point(1.0, d, 0)), p1)

    lc0 = clparam(ip.CubicCurve(pts=[Point(1-hl, mb, 0), Point(0.5, n0, 0)],
                                t=[0, 1], boundary=ip.TANGENT,
                                der=[.01 * Point(-sin(th), -cos(th)*N, 0),
                                     (.5-hl) * Point(-1, 0, 0)]), p1)
    lc1 = clparam(ip.CubicCurve(pts=[Point(1-hl, md, 0), Point(0.5, n2, 0)],
                                t=[0, 1], boundary=ip.TANGENT,
                                der=[.01 * Point(-sin(th), cos(th)*N, 0),
                                     (.5-hl) * Point(-1, 0, 0)]), p1)

    rc0 = clparam(ip.CubicCurve(pts=[Point(hl, mb, 0), Point(0.5, n0, 0)],
                                t=[0, 1], boundary=ip.TANGENT,
                                der=[.01 * Point(sin(th), -cos(th)*N, 0),
                                     (.5-hl) * Point(1, 0, 0)]), p1)
    rc1 = clparam(ip.CubicCurve(pts=[Point(hl, md, 0), Point(0.5, n2, 0)],
                                t=[0, 1], boundary=ip.TANGENT,
                                der=[.01 * Point(sin(th), cos(th)*N, 0),
                                     (.5-hl) * Point(1, 0, 0)]), p1)

    surfaces = [CoonsSurfacePatch([r2, mr1, flip(rt2), flip(rt3)]),
                CoonsSurfacePatch([l2, ml1, flip(lt2), flip(lt3)]),
                CoonsSurfacePatch([r1, mr1, flip(rt1), flip(rt0)]),
                CoonsSurfacePatch([l1, ml1, flip(lt1), flip(lt0)]),
                CoonsSurfacePatch([rt2, mr2, flip(t2), flip(rc1)]),
                CoonsSurfacePatch([lt2, ml2, flip(t2), flip(lc1)]),
                CoonsSurfacePatch([rt1, mr2, flip(t1), flip(rc0)]),
                CoonsSurfacePatch([lt1, ml2, flip(t1), flip(lc0)])]

    surfaces = map(hmSurf, surfaces)

    fix_surfs = [([r3, rt3, rc1, flip(t3)], 0.0, 0.5),
                 ([l3, lt3, lc1, flip(t3)], 0.5, 1.0),
                 ([r0, rt0, rc0, flip(t0)], 0.0, 0.5),
                 ([l0, lt0, lc0, flip(t0)], 0.5, 1.0)]
    for curves, ulow, uhigh in fix_surfs:
        presurf = CoonsSurfacePatch(map(hmCurve, curves))
        surfaces.append(ProjectSurfaceToSurface(presurf, height, ulow, uhigh))

    return surfaces


def WingTip(wv, params):

    def points(s):
        ku, kv, kw = s.GetKnots()
        ku = [ku[0], (ku[0]+ku[1])/2] + ku[1:-1] + [(ku[-2]+ku[-1])/2, ku[-1]]
        return [s.Evaluate(k, kv[0], kw[-1]) for k in ku]

    def derivatives(s):
        srf = s.GetFaces()[2]
        ku, kv = srf.GetKnots()
        ku = [ku[0], (ku[0]+ku[1])/2] + ku[1:-1] + [(ku[-2]+ku[-1])/2, ku[-1]]
        return [- (srf.EvaluateTangent(k, kv[-1])[0] %
                   srf.EvaluateNormal(k, kv[-1])).Normalize()
                for k in ku]

    def split_ub(s):
        btm = (s[0][:-2] + [s[0][-1]] +
               s[1][2:-2] + [s[1][-1]] +
               s[2][2:-2] + [s[2][-1]] +
               s[3][2:])
        top = (s[4][:-2] + [s[4][-1]] +
               s[5][2:-2] + [s[5][-1]] +
               s[6][2:-2] + [s[6][-1]] +
               s[7][2:])[::-1]
        return btm, top

    def refine_curve(s, N=100):
        k = s.GetKnots()

        kts = list(linspace(k[0], k[-1], N))
        kts = [kts[0], (kts[0]+kts[1])/2] + kts[1:-1] + [(kts[-2]+kts[-1])/2, kts[-1]]
        return [s.Evaluate(k) for k in kts]

    def heightMapSurface(height, surface):
        kus, kvs = surface.GetKnots()
        ekus = [kus[0], (kus[0]+kus[1])/2] + kus[1:-1] + [(kus[-2]+kus[-1])/2, kus[-1]]
        ekvs = [kvs[0], (kvs[0]+kvs[1])/2] + kvs[1:-1] + [(kvs[-2]+kvs[-1])/2, kvs[-1]]

        ipts = [surface.Evaluate(ku, kv) for ku in ekus for kv in ekvs]
        opts = [height.Evaluate(p[0], p[1]) for p in ipts]

        return ip.InterpolateSurface(opts, ekvs, ekus,
                                     [kvs[0]]*3 + kvs + [kvs[-1]]*3,
                                     [kus[0]]*3 + kus + [kus[-1]]*3).ReParametrize()

    def heightMapCurve(height, curve):
        okts = curve.GetKnots()
        ekts = [okts[0], (okts[0]+okts[1])/2] + okts[1:-1] + [(okts[-2]+okts[-1])/2, okts[-1]]

        ipts = [curve.Evaluate(k) for k in ekts]
        opts = [height.Evaluate(p[0], p[1]) for p in ipts]

        fkts = [okts[0]]*3 + okts + [okts[-1]]*3
        return ip.InterpolateCurve(opts, ekts, fkts)

    lines = map(points, wv)[::-1]; btm, top = split_ub(lines)
    derivs = map(derivatives, wv)[::-1]; btmd, topd = split_ub(derivs)

    topcurve, midp, lens = MkTopCurve(top, btm, topd, [(0.52, 0.2)])
    dist = ip.CubicCurve(pts=[Point(0,0,h) for h in lens], t=list(linspace(0, 1, len(lens))),
                         boundary=ip.NATURAL)

    pts = []

    for i, (t, b, td, bd, p) in enumerate(zip(top, btm, topd, btmd, midp)):
        if i == 0 or i == len(top) - 1:
            pts += [t]*102
        else:
            m = topcurve.Evaluate(p)
            md = (t - b).Normalize()

            h = m[2] - t[2]
            d = float(sqrt(h * abs(t - b)))

            bl = ip.CubicCurve(pts=[b, m], t=[0, 1], boundary=ip.TANGENT, der=[h*bd, d*md])
            br = ip.CubicCurve(pts=[t, m], t=[0, 1], boundary=ip.TANGENT, der=[h*td, -d*md])

            bl.AppendCurve(br.FlipParametrization(), continuity=1, reparam=True)
            pts += refine_curve(bl.FlipParametrization(), 100)


    vkts = map(float, range(len(top)-2))
    vkts = [vkts[0], (vkts[0]+vkts[1])/2] + vkts[1:-1] + [(vkts[-2]+vkts[-1])/2, vkts[-1]]
    vkts_ = [vkts[0]]*4 + vkts[2:-2] + [vkts[-1]]*4
    ekts = list(np.linspace(0, 1, 100))
    ekts = [ekts[0], (ekts[0]+ekts[1])/2] + ekts[1:-1] + [(ekts[-2]+ekts[-1])/2, ekts[-1]]
    ekts_ = [ekts[0]]*4 + ekts[2:-2] + [ekts[-1]]*4

    height = ip.InterpolateSurface(pts, ekts, vkts, ekts_, vkts_)

    N1, N2 = map(lambda l: len(l) - 3, lines[:2])
    distfun = lambda d: (N1+N2) * 2 * dist.GetParameterAtPoint(Point(0,0,d))[0]

    surfaces = MkSurfaces(N1, N2, params, distfun,
                          lambda s: heightMapSurface(height, s),
                          lambda c: heightMapCurve(height, c),
                          height)

    WriteG2('out/test.g2', surfaces + wv)

    return surfaces


def OMeshTip(z, R, cs, tip):

    center = Point(0, 0, z)

    def pillar(p, ang, N):
        axis = (p - center) % center

        return CSeg(center, p, ang, axis, N).ReParametrize()

    def summit(p, N):
        axis = (p.Evaluate(0.0) - center) % center
        return CSeg(center, p.Evaluate(1.0), pi/4, axis, N).ReParametrize()

    def beam(o, p1, p2, sign, N):
        start = p1.Evaluate(1.0)
        end = p2.Evaluate(1.0)
        axis = sign * ((start - center) % (end - center))
        angle = arccos((start - center).Normalize() * (end - center).Normalize())
        return CSeg(center, start, sign * angle, axis, N).ReParametrize()

    def icirc(start, end, N, space='lin'):
        axis = ((start - center) % (end - center)).Normalize()
        angle = arccos((end - center).Normalize() * (start - center).Normalize())
        return CSeg(center, start, angle, axis, N, space=space).ReParametrize()

    def mkWall(outer, beam, N):
        icircs = [icirc(outer.Evaluate(k), beam.Evaluate(k), N, space='lin')
                  for k in outer.GetKnots()]
        pts = list(chain(*[[ic.Evaluate(k) for k in ic.GetKnots()] for ic in icircs]))
        return InterpolateSurface(pts,
                                  list(linspace(0, 1, N+1)),
                                  list(linspace(0, 1, len(outer.GetKnots()))))

    def mkCeiling(o1, o2, beam, summit):
        icircs = [icirc(beam.Evaluate(k), summit.Evaluate(k), len(o2.GetKnots()) - 1)
                  for k in o1.GetKnots()]
        pts = list(chain(*[[ic.Evaluate(k) for k in ic.GetKnots()] for ic in icircs]))
        return InterpolateSurface(pts,
                                  list(linspace(0, 1, len(o2.GetKnots()))),
                                  list(linspace(0, 1, len(o1.GetKnots()))))

    outer = [v.GetEdges()[2].ReParametrize() for v in cs]
    rays = [v.GetEdges()[3].ReParametrize() for v in cs]
    pillars = [pillar(r.Evaluate(1.0), a, 20) for r, a in zip(rays, cycle([pi/4, arccos(sqrt(2.0/3.0))]))]

    beams = [beam(outer[0], pillars[0], pillars[1], -1, 20), beam(outer[1], pillars[2], pillars[1], 1, 20),
             beam(outer[2], pillars[2], pillars[3], -1, 20), beam(outer[3], pillars[4], pillars[3], 1, 20),
             beam(outer[4], pillars[4], pillars[5], -1, 20), beam(outer[5], pillars[6], pillars[5], 1, 20),
             beam(outer[6], pillars[6], pillars[7], -1, 20), beam(outer[7], pillars[0], pillars[7], 1, 20)]

    summits = [summit(p, 20).FlipParametrization() for p in pillars[::2]]

    for o in outer[1::2]:
        o.FlipParametrization()

    N = len(tip[0].GetKnots()[0]) - 1

    walls = [mkWall(o, b, N) for o, b in zip(outer, beams)]

    ceilings = [mkCeiling(outer[0], outer[1], beams[0], summits[1]),
                mkCeiling(outer[2], outer[3], beams[2], summits[2]),
                mkCeiling(outer[4], outer[5], beams[4], summits[3]),
                mkCeiling(outer[6], outer[7], beams[6], summits[0])]

    return walls + ceilings


def TipTFISurfaces(cs, tip, sphere, params):

    def DoLine(source, source_n, target, target_n, r=None, fknots=None):
        crv = InterpolateCurve([source, target], [0, abs(source-target)], source_n, target_n)

        clen = abs(source - target)
        coeff = (1.0 - pow(r, params.Nrad)) / (1.0 - r)
        knots = GradedSpace(0.0, clen / coeff, r, params.Nrad) + [clen]
        knots = [knots[0], 0.5 * (knots[0] + knots[1])] + knots[1:]
        knots = [k/clen*abs(source-target) for k in knots]
        pts = [crv.Evaluate(k) for k in knots]

        return InterpolateCurve(pts, fknots)


    def MkWallPatch(tips, spheres, cs, line, flip=False):
        (ku1, kv1), (ku2, kv2), (ku3, kv3), (ku4, kv4) = [t.GetKnots() for t in tips]

        ipts = (GetCurvePoints(tips[0].GetEdges()[2 if flip else 3]) +
                GetCurvePoints(tips[2].GetEdges()[2 if flip else 1])[1:])
        inner = InterpolateCurve(ipts, range(len(ipts)))

        opts = GetCurvePoints(spheres[0].GetEdges()[0]) + GetCurvePoints(spheres[1].GetEdges()[0])[1:]
        outer = InterpolateCurve(opts, range(len(opts)))

        if flip:
            normals = [(tips[0].EvaluateNormal(k1, kv1[-1]) -
                        tips[1].EvaluateNormal(k2, kv2[-1])).Normalize()
                       for k1, k2 in zip(ku1, ku2)]
            split = len(normals)
            normals += [(tips[2].EvaluateNormal(k3, kv3[-1]) -
                         tips[3].EvaluateNormal(k4, kv4[-1])).Normalize()
                        for k3, k4 in zip(ku3, ku4)]
        else:
            normals = [(tips[0].EvaluateNormal(ku1[0], k1) -
                        tips[1].EvaluateNormal(ku2[0], k2)).Normalize()
                       for k1, k2 in zip(kv1, kv2)]
            split = len(normals)
            normals += [(tips[3].EvaluateNormal(ku4[-1], k4) -
                         tips[2].EvaluateNormal(ku3[-1], k3)).Normalize()
                        for k3, k4 in zip(kv3, kv4)[1:]]

        srf = TFI.OrthogonalSurface([inner, outer, cs, line], normals,
                                    ranges=[(0, len(normals))], fac_blend=0.9,
                                    bnd_layer=params.bndlayer)
        kus, kvs = srf.GetKnots()
        return [srf.GetSubSurf((kus[0], kvs[0]), (kus[split-1], kvs[-1])),
                srf.GetSubSurf((kus[split-1], kvs[0]), (kus[-1], kvs[-1]))]


    def MkSmallWallPatch(t1, t2, pillar, cs, line, flip=False):
        (ku1, kv1), (ku2, kv2) = [t.GetKnots() for t in [t1,t2]]

        normals = [(t1.EvaluateNormal(ku1[-1], k1) - t2.EvaluateNormal(k2, kv2[0])).Normalize()
                   * (1 if flip else -1)
                   for k1, k2 in zip(kv1, ku2)]

        return TFI.OrthogonalSurface([t1.GetEdges()[1], pillar, cs, line], normals,
                                     ranges=[(0, len(normals))], fac_blend=0.9,
                                     bnd_layer=params.bndlayer)


    def MkFloorPatches(t1, t2, t3, s1, s2, w1, w2, line, flip=False):
        (ku1, kv1), (ku2, kv2), (ku3, kv3) = [t.GetKnots() for t in [t1,t2,t3]]

        normals = [(t3.EvaluateNormal(k3, kv3[0]) - t1.EvaluateNormal(k1, kv1[-1])).Normalize()
                   * (-1 if flip else 1) for k3, k1 in zip(ku3[::-1], ku1)]
        srf1 = TFI.OrthogonalSurface([t1.GetEdges()[2], s1.GetEdges()[1], w1, line], normals,
                                     ranges=[(0, len(normals))], fac_blend=0.85,
                                     bnd_layer=params.bndlayer)

        normals = [(t2.EvaluateNormal(ku2[-1], k2) + t3.EvaluateNormal(ku3[0], k3)).Normalize()
                   * (-1 if flip else 1) for k2, k3 in zip(kv2, kv3)]
        srf2 = TFI.OrthogonalSurface([t2.GetEdges()[1], s2.GetEdges()[1].FlipParametrization(), line, w2],
                                     normals, ranges=[(0, len(normals))], fac_blend=0.85,
                                     bnd_layer=params.bndlayer)

        return [srf1, srf2]


    def CenterLine(t1, t2, t3, s1, s2, s3, c, w1, w2, flip=False):
        def Line(tp, sp, l, r, ev):
            c1 = InterpolateCurve(tp, range(len(tp)), order=2)
            c2 = InterpolateCurve(sp, range(len(sp)), order=2)
            return TFI.LinearSurface([c1, c2, l, r], eval_xi=[ev], interpolate=False)

        if flip:
            s3 = s3.Clone().SwapParametrization()
            s3.FlipParametrization(0).FlipParametrization(1)

        (ku1, kv1), (ku2, kv2), (ku3, kv3) = [t.GetKnots() for t in [t1,t2,t3]]

        tp = GetCurvePoints(t1.GetEdges()[1]) + GetCurvePoints(t3.GetEdges()[0])[1:]
        sp = GetCurvePoints(s1.GetEdges()[2]) + GetCurvePoints(s1.GetEdges()[1])[-1::-1]
        pts1 = Line(tp, sp, c, w1, len(ku2) - 1)

        tp = GetCurvePoints(t2.GetEdges()[0]) + GetCurvePoints(t2.GetEdges()[1])[1:]
        sp = GetCurvePoints(s2.GetEdges()[2]) + GetCurvePoints(s3.GetEdges()[2])[1:]
        pts2 = Line(tp, sp, c, w2, len(ku2) - 1)

        tp = GetCurvePoints(t1.GetEdges()[2]) + GetCurvePoints(t2.GetEdges()[1])[1:]
        sp = GetCurvePoints(s3.GetEdges()[3]) + GetCurvePoints(s3.GetEdges()[2])[1:]
        pts3 = Line(tp, sp, w1, w2, len(ku3) - 1)

        pts = [(p1 + p2 + p3) / 3 for p1, p2, p3 in zip(pts1, pts2, pts3)]
        ikts = CurveLengthParametrization(pts, True)
        crv = InterpolateCurve(pts, range(len(pts)))

        N = 8*params.bndlayer

        nrm = (t1.EvaluateNormal(ku1[-1], kv1[-1]) -
               t2.EvaluateNormal(ku2[-1], kv2[0]) -
               t3.EvaluateNormal(ku3[0], kv3[0])).Normalize() * (1 if flip else -1)
        tnrm = crv.EvaluateTangent(N).Normalize()
        pt, tpt = crv.Evaluate(0.0), crv.Evaluate(N)
        dist = abs(pt - tpt)

        tcrv = InterpolateCurve([pt, tpt], [0.0, 1.0], 0.5*dist*nrm, dist*tnrm, order=4)
        UniformCurve(tcrv, n=100)
        tpts = GetCurvePoints(tcrv) + pts[N+1:]
        tcrv = InterpolateCurve(tpts, CurveLengthParametrization(tpts, True))

        pts = [tcrv.Evaluate(k) for k in ikts]

        return InterpolateCurve(pts, range(len(pts)))


    ku, kv = tip[7].GetKnots()
    source = tip[7].Evaluate(ku[-1], kv[-1])
    source_n = tip[7].EvaluateNormal(ku[-1], kv[-1])
    ku, kv = tip[6].GetKnots()
    source_n = (source - tip[6].EvaluateNormal(ku[-1], kv[-1])).Normalize()

    ku, kv = sphere[8].GetKnots()
    target = sphere[8].Evaluate(ku[-1], kv[0])
    target_n = sphere[8].EvaluateNormal(ku[-1], kv[0])

    line = DoLine(source, source_n, target, target_n, 1.0/0.9, cs[0].GetKnots()[1])

    walls =  MkWallPatch([tip[k] for k in [10,11,6,7]], [sphere[k] for k in [0, 8]],
                         cs[0].GetEdges()[3], line)
    walls += MkWallPatch([tip[k] for k in [3,1,7,5]], [sphere[k] for k in [2, 9]],
                         cs[2].GetEdges()[3], line, flip=True)
    walls += MkWallPatch([tip[k] for k in [9,8,5,4]], [sphere[k] for k in [4, 10]],
                         cs[4].GetEdges()[3], line)
    walls += MkWallPatch([tip[k] for k in [0,2,4,6]], [sphere[k] for k in [6, 11]],
                         cs[6].GetEdges()[3], line, flip=True)

    centers = [CenterLine(tip[11], tip[3], tip[7], sphere[0], sphere[1], sphere[8],
                          cs[1].GetEdges()[3], walls[0].GetEdges()[1], walls[2].GetEdges()[1]),
               CenterLine(tip[9], tip[1], tip[5], sphere[3], sphere[2], sphere[9],
                          cs[3].GetEdges()[3], walls[4].GetEdges()[1], walls[2].GetEdges()[1], True),
               CenterLine(tip[8], tip[0], tip[4], sphere[4], sphere[5], sphere[10],
                          cs[5].GetEdges()[3], walls[4].GetEdges()[1], walls[6].GetEdges()[1]),
               CenterLine(tip[10], tip[2], tip[6], sphere[7], sphere[6], sphere[11],
                          cs[7].GetEdges()[3], walls[0].GetEdges()[1], walls[6].GetEdges()[1], True)]

    smallwalls = [MkSmallWallPatch(tip[11], tip[3], sphere[0].GetEdges()[2],
                                   cs[1].GetEdges()[3], centers[0]),
                  MkSmallWallPatch(tip[9], tip[1], sphere[2].GetEdges()[2],
                                   cs[3].GetEdges()[3], centers[1], True),
                  MkSmallWallPatch(tip[8], tip[0], sphere[4].GetEdges()[2],
                                   cs[5].GetEdges()[3], centers[2]),
                  MkSmallWallPatch(tip[10], tip[2], sphere[6].GetEdges()[2],
                                   cs[7].GetEdges()[3], centers[3], True)]

    floors =  MkFloorPatches(tip[11], tip[3], tip[7], sphere[0], sphere[1],
                             walls[0].GetEdges()[1], walls[2].GetEdges()[1], centers[0])
    floors += MkFloorPatches(tip[9], tip[1], tip[5], sphere[3], sphere[2],
                             walls[4].GetEdges()[1], walls[2].GetEdges()[1], centers[1], True)
    floors += MkFloorPatches(tip[8], tip[0], tip[4], sphere[4], sphere[5],
                             walls[4].GetEdges()[1], walls[6].GetEdges()[1], centers[2])
    floors += MkFloorPatches(tip[10], tip[2], tip[6], sphere[7], sphere[6],
                             walls[0].GetEdges()[1], walls[6].GetEdges()[1], centers[3], True)

    return walls + smallwalls + floors



def TipCircleVolumes(cs, tip, sphere, tfi):

    def ImposeGrading(srf):
        srf.ReParametrize()
        kus, kvs = srf.GetKnots()

        curves = [InterpolateCurve([srf.Evaluate(ku, kv) for ku in kus], kus) for kv in kvs]

        N = len(kus)
        ds = 0.008
        fac = ComputeFactor(1.0, ds, N-1)
        graded = GradedSpace(0.0, ds, fac, N)
        graded[-1] = 1.0

        rs = list(linspace(0, 1, len(kvs)))
        rs = [max(arctan(1.0*(r - 0.5))/pi, 0.0) for r in rs]
        rs = [r/rs[-1] for r in rs]

        def GradeCurve(c, r):
            space = [r*g + (1.0-r)*s for g, s in zip(graded, kus)]
            return InterpolateCurve([c.Evaluate(s) for s in space], kus, order=2)

        curves = [GradeCurve(c, r) for c, r in zip(curves, rs)]

        return LoftCurves(curves, order=2)


    def TFIVolumeActual(s1, s2, s3, s4, s5, s6, flipu=False, flipv=False, grad=True):
        for s in [s1, s2, s3, s4, s5, s6]:
            s.ReParametrize()

        if flipu:
            s1 = s1.Clone().FlipParametrization(0)
            s2 = s2.Clone().FlipParametrization(0)

        if flipv:
            s5 = s5.Clone().FlipParametrization(0)
            s6 = s6.Clone().FlipParametrization(0)

        kus, kvs = s1.GetKnots()
        _, kws = s3.GetKnots()

        rws = CurveLengthParametrization([s3.Evaluate(kus[len(kus)/2], w) for w in kws], normalize=True)

        surfs = []

        # Yes, this is slow.
        for w, rw in zip(kws, rws):
            pts = [(1-rw)*s1.Evaluate(u,v) + rw*s2.Evaluate(u,v) +
                   (1-v)*s3.Evaluate(u,w) + v*s4.Evaluate(u,w) +
                   (1-u)*s5.Evaluate(v,w) + u*s6.Evaluate(v,w) -
                   ((1-u)*(1-v)*s3.Evaluate(0,w) + u*v*s4.Evaluate(1,w) +
                    u*(1-v)*s3.Evaluate(1,w) + (1-u)*v*s4.Evaluate(0,w)) -
                   ((1-v)*(1-rw)*s1.Evaluate(u,0) + v*rw*s2.Evaluate(u,1) +
                    v*(1-rw)*s1.Evaluate(u,1) + (1-v)*rw*s2.Evaluate(u,0)) -
                   ((1-u)*(1-rw)*s1.Evaluate(0,v) + u*rw*s2.Evaluate(1,v) +
                    u*(1-rw)*s1.Evaluate(1,v) + (1-u)*rw*s2.Evaluate(0,v)) +
                   (1-u)*(1-v)*(1-rw)*s1.Evaluate(0,0) + u*v*rw*s2.Evaluate(1,1) +
                   u*(1-v)*(1-rw)*s1.Evaluate(1,0) + u*v*(1-rw)*s1.Evaluate(1,1) +
                   (1-u)*(1-v)*rw*s2.Evaluate(0,0) + (1-u)*v*rw*s2.Evaluate(0,1) +
                   (1-u)*v*(1-rw)*s1.Evaluate(0,1) + u*(1-v)*rw*s2.Evaluate(1,0)
                   for v in kvs for u in kus]
            surfs.append(InterpolateSurface(pts, kus, kvs))

        return LoftSurfaces(surfs, kws, order=4)

    tip = [t.Clone() for t in tip]
    sphere = [s.Clone() for s in sphere]
    tfi = [t.Clone() for t in tfi]

    sphere[0].SwapParametrization()
    sphere[1].SwapParametrization().FlipParametrization(0)
    sphere[2].SwapParametrization()
    sphere[3].SwapParametrization().FlipParametrization(0)
    sphere[4].SwapParametrization()
    sphere[5].SwapParametrization().FlipParametrization(0)
    sphere[6].SwapParametrization()
    sphere[7].SwapParametrization().FlipParametrization(0)
    sphere[8].FlipParametrization(1)
    sphere[9].SwapParametrization().FlipParametrization(0)
    sphere[10].FlipParametrization(1)
    sphere[11].SwapParametrization().FlipParametrization(0)

    tip[3].SwapParametrization()
    tip[1].SwapParametrization().FlipParametrization(0)
    tip[9].FlipParametrization(0)
    tip[0].SwapParametrization()
    tip[2].SwapParametrization().FlipParametrization(0)
    tip[10].FlipParametrization(0)
    tip[7].SwapParametrization()
    tip[5].SwapParametrization()
    tip[4].SwapParametrization()
    tip[6].SwapParametrization()

    tfi[15].FlipParametrization(0)
    tfi[14].FlipParametrization(0)
    tfi[19].FlipParametrization(0)
    tfi[18].FlipParametrization(0)

    vols = [TFIVolumeActual(tip[11], sphere[0], cs[0], tfi[12], tfi[0], tfi[8]),
            TFIVolumeActual(tip[3], sphere[1], cs[1], tfi[13], tfi[8], tfi[2]),
            TFIVolumeActual(tip[1], sphere[2], cs[2], tfi[15], tfi[2], tfi[9]),
            TFIVolumeActual(tip[9], sphere[3], cs[3], tfi[14], tfi[9], tfi[4]),
            TFIVolumeActual(tip[8], sphere[4], cs[4], tfi[16], tfi[4], tfi[10]),
            TFIVolumeActual(tip[0], sphere[5], cs[5], tfi[17], tfi[10], tfi[6]),
            TFIVolumeActual(tip[2], sphere[6], cs[6], tfi[19], tfi[6], tfi[11]),
            TFIVolumeActual(tip[10], sphere[7], cs[7], tfi[18], tfi[11], tfi[0])]

    tfi[12].FlipParametrization(0)
    tfi[15].FlipParametrization(0)
    tfi[16].FlipParametrization(0)
    tfi[19].FlipParametrization(0)

    vols += [TFIVolumeActual(tip[7], sphere[8], tfi[13], tfi[1], tfi[12], tfi[3]),
             TFIVolumeActual(tip[5], sphere[9], tfi[15], tfi[5], tfi[14], tfi[3]),
             TFIVolumeActual(tip[4], sphere[10], tfi[17], tfi[5], tfi[16], tfi[7]),
             TFIVolumeActual(tip[6], sphere[11], tfi[19], tfi[1], tfi[18], tfi[7])]

    return vols



def OuterVolumes(os, tipvols, params):

    def icirc(pa, pb, ikts, okts):
        h = pb[2] - pa[2]
        proj = Point(pb[0], pb[1], pa[2])
        d = abs(proj - pa)
        L = (d**2 + h**2) / (2*d)
        center = pa + (proj - pa).Normalize() * L

        cseg = CSeg(center, pa, arccos((pa-center).Normalize() * (pb-center).Normalize()),
                    ((pa - center) % (pb - center)).Normalize(), len(ikts) - 1).ReParametrize()
        return InterpolateCurve([cseg.Evaluate(k) for k in ikts], okts, order=2)

    def MkLowerWall(inner, btm):
        h = 2 * params.R
        outer = InterpolateCurve([btm.Evaluate(btm.GetKnots()[-1]),
                                  btm.Evaluate(btm.GetKnots()[-1]) + Point(0,0,h)], [0.0, 1.0], order=2)
        ikts = UniformParametrization(GetCurvePoints(inner), normalize=True)
        outer = InterpolateCurve([outer.Evaluate(k) for k in ikts], inner.GetKnots(), order=2)

        top = InterpolateCurve([inner.Evaluate(1.0), outer.Evaluate(1.0)], [0.0, 1.0], order=2)
        ikts = CurveLengthParametrization(GetCurvePoints(btm), normalize=True)
        top = InterpolateCurve([top.Evaluate(k) for k in ikts], btm.GetKnots())

        ikts = UniformParametrization(GetCurvePoints(inner), normalize=True)
        okts = inner.GetKnots()
        edges = [icirc(btm.Evaluate(k), top.Evaluate(k), ikts, okts) for k in top.GetKnots()[1:-1]]
        edges = [inner] + edges + [outer]

        return LoftCurves(edges, btm.GetKnots(), order=2)

    def MkLowerVolume(inner, osrf):
        bu, bv = osrf.GetKnots()
        iu, iv, iw = inner.GetKnots()
        iw = iw[-1]

        walls = []
        for bp, ip in zip(bu, iu):
            btm = InterpolateCurve([osrf.Evaluate(bp, k) for k in bv], bv, order=2)
            inr = InterpolateCurve([inner.Evaluate(ip, k, iw) for k in iv], iv, order=2)
            walls.append(MkLowerWall(inr, btm))

        return LoftSurfaces(walls, iu, order=2)

    def MkTopVolume(inner, left, right, ikts, flip=False):
        top = Point(0, 0, params.wingdef[-1].z + 2 * params.R)

        e1 = left.GetFaces()[1].GetEdges()[1].ReParametrize()
        e2 = right.GetFaces()[1].GetEdges()[1].ReParametrize()

        if flip:
            e1.FlipParametrization()
            e2.FlipParametrization()

        kts1 = list(linspace(0, 1, len(e1.GetKnots())))
        kts2 = list(linspace(0, 1, len(e2.GetKnots())))
        e3 = InterpolateCurve([top, e2.Evaluate(1.0)], [0.0, 1.0], order=2)
        e3 = InterpolateCurve([e3.Evaluate(k) for k in kts1], e1.GetKnots(), order=2).FlipParametrization()
        e4 = InterpolateCurve([e1.Evaluate(0.0), top], [0.0, 1.0], order=2)
        e4 = InterpolateCurve([e4.Evaluate(k) for k in kts2], e2.GetKnots(), order=2).FlipParametrization()

        ceiling = CoonsSurfacePatch([e1, e2, e3, e4])
        ceiling.SwapParametrization()

        inner = inner.Clone().FlipParametrization(1)

        rus, rvs, rws = right.GetKnots()
        lus, lvs, lws = left.GetKnots()
        ius, ivs, iws = inner.GetKnots(); iw = iws[-1]
        cus, cvs = ceiling.GetKnots()

        surfaces = []
        for iu, cu in zip(ius, cus):
            edges = []
            for iv, cv in zip(ivs, cvs):
                ip = inner.Evaluate(iu, iv, iw)
                cp = ceiling.Evaluate(cu, cv)
                edges.append(InterpolateCurve([(1-k)*ip + k*cp for k in ikts], ikts))
            surfaces.append(LoftCurves(edges, e1.GetKnots(), order=2))

        volume = LoftSurfaces(surfaces, e2.GetKnots(), order=2)

        return volume

    lowerVols = [MkLowerVolume(i, o) for i, o in zip(tipvols[:8], os)]

    lu, lv, lw = lowerVols[0].GetKnots()
    ikts = CurveLengthParametrization([lowerVols[0].Evaluate(lu[0], k, lw[0]) for k in lv], normalize=True)

    topVols = [MkTopVolume(tipvols[8], lowerVols[0], lowerVols[1], ikts),
               MkTopVolume(tipvols[9], lowerVols[3], lowerVols[2], ikts, flip=True),
               MkTopVolume(tipvols[10], lowerVols[4], lowerVols[5], ikts),
               MkTopVolume(tipvols[11], lowerVols[7], lowerVols[6], ikts, flip=True)]

    return lowerVols + topVols
