# Definition of sections for NREL 5 MW turbine
# References
#
#  "Definition of a 5-MW Reference Wind Turbine
#  for Offshore System Development". J. Jonkman,
#  S. Butterfield, W.Musial and G, Scott. Technical
#  Report NREL/TP-500-38060, February 2009.
#
#  "Development of a Scale Model for Testing of Offshore
#  Floating Wind Turbine systems". H. R. Martin, MSc thesis,
#  The University of Maine, December, 2011.

from GoTools import *
# from GoTools.CurveFactory import *
# from GoTools.SurfaceFactory import *
from GoTools.VolumeFactory import *
from GoTools.Preprocess import *
from GeoUtils.CurveUtils import *
from GeoUtils.Elementary import *
from GeoUtils.IO import *
from GeoUtils.Refinement import *
import GeoUtils.Interpolate as ip
import GeoUtils.TFI as tfi

from Rounding import *

from collections import namedtuple
from itertools import repeat, izip_longest, chain
from math import *
import numpy as np
from operator import itemgetter, attrgetter
import shutil
import os
import csv
import sys
import xml.etree.ElementTree as ET


edgemap = {0: [9,7,11,5],
           1: [10,8,12,6],
           2: [1,6,2,5],
           3: [3,8,4,7],
           4: [1,10,3,9],
           5: [2,12,4,11]}

def GetEdge(vol, i):
    for f_idx, e_idxs in edgemap.iteritems():
        if i in e_idxs:
            return vol.GetFaces()[f_idx].GetEdges()[e_idxs.index(i)]


def PrepareParams(params):
    for arg in sys.argv:
        for key in params.keys():
            if arg in [key, 'no' + key] and type(params[key]) == bool:
                params[key] = arg == key
            elif arg.startswith(key + '='):
                arg = arg[len(key)+1:]
                if type(params[key]) in [float, int]:
                    params[key] = type(params[key])(arg)
                else:
                    params[key] = arg

    def fix(d):
        for k in ['z', 'theta', 'chord', 'ac', 'ao']:
            d[k] = float(d[k])
        return d

    wingdef_tree = ET.parse(params['wingfile'])
    params['wingdef'] = [namedtuple('Section', s.attrib.keys())(**fix(s.attrib))
                         for s in wingdef_tree.getroot()]

    params['teC'] = params['teC_fac'] * params['te']
    params['computed'] = {}

    return namedtuple('Params', params.keys())(**params)


def PrepareOutput(params):
    shutil.rmtree('out', ignore_errors=True)

    folders = ['out']
    if params.debug:
        folders += map(lambda s: 'out/'+s, ['ndsections', 'sections', 'wingsecs',
                                            'circlesecs', 'tfi',      'crossecs',
                                            'inner',      'outer',    'sections2d',
                                            'bnds',       'groups'])
    for folder in folders:
        try: os.makedirs(folder)
        except OSError: pass


def PrintPatchSizes(revols):
    ndofs = {}
    for i, v in enumerate(revols):
        ku, kv, kw = v.GetKnots()
        n = len(ku) * len(kv) * len(kw)
        if n not in ndofs:
            ndofs[n] = []
        ndofs[n].append(i+1)

    def ranges(lst):
        split = []
        for i, (prv, nxt) in enumerate(zip(lst[:-1], lst[1:])):
            if nxt - prv > 1:
                split.append(i+1)
        rngs = []
        split = [None] + split + [None]
        for prv, nxt in zip(split[:-1], split[1:]):
            lstpart = lst[slice(prv, nxt)]
            rngs.append((lstpart[0], lstpart[-1]))
        return rngs

    print ''

    for n in sorted(ndofs, reverse=True):
        s = []
        for a,b in ranges(ndofs[n]):
            s.append(('%i-%i' % (a,b)) if a != b else ('%i' % a))

        title = '%i DoFs:' % n
        indent = len(title)

        print title,

        n = 0
        for i, p in enumerate(s):
            if n + len(p) + 2 > 60:
                print ''
                print ' '*indent,
                n = 0
            print (p + ',') if i < len(s)-1 else p,
            n += len(p) + 2

        print ''

    print 'Total: ~%i DoFs' % sum(n * len(patches) for n, patches in ndofs.iteritems())


def WriteXML(patchgroups, face_bnds, edge_bnds, reorder, params):
    old_to_new = lambda i: reorder.index(i)
    new_to_old = lambda i: reorder[i]

    ends = cumsum(patchgroups)
    starts = hstack(([1], ends[:-1] + 1))

    pg = ['    <part proc="%i" lower="%i" upper="%i" />' % tup
          for tup in zip(range(len(starts)), starts, ends)]

    ts_face = [['    <set name="%s" type="face">' % name]
             + ['      <item patch="%i">%i</item>' % (old_to_new(pid-1)+1, face) for pid, face in boundary]
             + ['    </set>']
               for name, boundary in face_bnds.iteritems()]
    ts_edge = [['    <set name="%s" type="edge">' % name]
             + ['      <item patch="%i">%i</item>' % (old_to_new(pid-1)+1, edge) for pid, edge in boundary]
             + ['    </set>']
               for name, boundary in edge_bnds.iteritems()]
    ts = list(chain(*(ts_face + ts_edge)))

    with open('out/%s.xinp' % params.out, 'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n\n')
        f.write('<geometry>\n')
        f.write('  <partitioning procs="%i">\n' % len(pg))
        for s in pg:
            f.write(s + '\n')
        f.write('  </partitioning>\n')
        f.write('  <patchfile>%s.g2</patchfile>\n'%params.out)
        f.write('  <nodefile>%s.hdf5</nodefile>\n'%params.out)
        f.write('  <topologysets>\n')
        for s in ts:
            f.write(s + '\n')
        f.write('  </topologysets>\n')
        f.write('</geometry>')


def CurveLength(pts):
    return sum([abs(pn-pp) for pp, pn in zip(pts[:-1], pts[1:])])


# Compute geometric refinement factor from curvelength, starting cellsize and number of cells
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


def F(s, dx1, dx2, n1, n2, r10, r20):
    r2 = r20
    n1i = 1.0 / n1
    coeff1 = dx1 / dx2
    r2n2 = pow(r2, n2)
    r1 = pow(coeff1*r2n2, n1i)
    dx = dx2 / r2n2

    return s - dx * (
        (1.0-pow(r1,n1+1)) / (1.0-r1) + (1.0-pow(r2,n2+1)) / (1.0-r2)
    )


def dFdr2(dx1, dx2, n1, n2, r10, r20):
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


def GradingTwoSided(s, dx1, dx2, n1, n2, r10, r20, tol, maxit):
    eps = 10.0 * tol
    r2 = r20
    r1 = pow(dx1 / dx2 * pow(r2,n2), 1.0/n1)

    it = 0
    while eps > tol and it <= maxit:
        P = F(s,dx1,dx2,n1,n2,r1,r2)
        dP = dFdr2(dx1,dx2,n1,n2,r1,r2)
        dr2 = -P / dP
        r2 = r2 + dr2
        r1 = pow(dx1 / dx2 * pow(r2,n2), 1.0/n1)
        eps = abs(P)
        it += 1

    if it == maxit:
        r1 = 0.0
        r2 = 0.0

    dx = dx2 / pow(r2,n2)

    return r1, r2, dx


def LoadWingSectionFtTXT(filen):
    # Load csv file
    data = csv.reader(open(filen,'rb'))
    rows = [[float(x) for x in r[0].split()] for r in data]
    return [map(itemgetter(i), rows) for i in range(4)] # xc, tc, f, dfdx


def ResampleWing(surface, params):
    knots1, knots2 = surface.GetKnots()

    critical = params.wingdef[params.addpt].z

    ylen = critical - knots2[0]
    ry = ComputeFactor(ylen, params.dzJoin, params.NlenBase)
    knoty = GradedSpace(critical, -params.dzJoin, ry, params.NlenBase + 1)[::-1]

    ylen = knots2[-1] - knoty[-1]
    ry1, ry2, _ = GradingTwoSided(ylen, params.dzJoin, params.dzTip,
                                  params.Nlen, params.Nlen, 0.9, 0.9, 1e-12, 200)
    knoty += GradedSpace(knoty[-1], params.dzJoin, 1.0/ry1, params.Nlen+2)[1:]
    knoty += GradedSpace(knots2[-1], -params.dzTip, 1.0/ry2, params.Nlen+1)[::-1]

    # Resample in Ny sections along the wing
    sections = []
    for eta in knoty:
        pts = [surface.Evaluate(k,eta) for k in knots1]
        sections.append(ip.CubicCurve(pts=pts, t=knots1, boundary=ip.PERIODIC))

    return sections


def TrailingEdgeModificationFt(sec, te):
    xc, tc, f, dfdx = LoadWingSectionFtTXT(sec)

    # Number of points in section
    np = len(xc)

    # X, Y coordinates of section
    # f and dfdx gives a function and its derivative of xc, which defines the center line.
    # tc gives the thickness normal to the center line at each x-point.
    x2D_u = [xci - 0.5 * tci * sin(atan(dfdxi)) for xci, tci, dfdxi in zip(xc, tc, dfdx)]
    y2D_u = [fi  + 0.5 * tci * cos(atan(dfdxi)) for fi,  tci, dfdxi in zip(f,  tc, dfdx)]
    x2D_l = [xci + 0.5 * tci * sin(atan(dfdxi)) for xci, tci, dfdxi in zip(xc, tc, dfdx)]
    y2D_l = [fi  - 0.5 * tci * cos(atan(dfdxi)) for fi,  tci, dfdxi in zip(f,  tc, dfdx)]

    # Midpoints
    mids = [Point((xu+xl)/2, (yu+yl)/2, 0) for xu, xl, yu, yl in zip(x2D_u, x2D_l, y2D_u, y2D_l)]
    params = CurveLengthParametrization(mids, normalize=True)

    # The data files already include a gap at the trailing edge.
    # Move the profiles slightly so that the trailing edge gap is equal to te.
    diff = (y2D_u[-1] - y2D_l[-1] - te) / 2
    y2D_u = [y - p*diff for y, p in zip(y2D_u, params)]
    y2D_l = [y + p*diff for y, p in zip(y2D_l, params)]

    x2D_nd = x2D_u[::-1] + x2D_l[1:]
    y2D_nd = y2D_u[::-1] + y2D_l[1:]

    return x2D_nd, y2D_nd


def TrailingEdgeCurve(wingcurve):
    # Create rounded trailing edge curve
    knots = wingcurve.GetKnots()

    p1 = wingcurve.Evaluate(knots[0])
    t1 = wingcurve.EvaluateTangent(knots[0]).Normalize()

    pN = wingcurve.Evaluate(knots[-1])
    tN = wingcurve.EvaluateTangent(knots[-1]).Normalize()
    pm = 0.5 * (p1 + pN)

    ttelen = abs(p1 - pN)
    n = (tN - t1).Normalize()
    pte = pm + 0.5 * ttelen * n

    crv = ip.CubicCurve(pts=[pN, pte, p1], boundary=ip.TANGENT, der=[tN, t1])
    UniformCurve(crv, n=5)

    return GetCurvePoints(crv), crv


def RotateCoordinates(xcoords, ycoords, x0, theta):
    np = len(xcoords)
    xc, yc = x0[0], x0[1]

    # Convert from degrees to radians
    alpha = pi * theta / 180.0

    xrot = [xc + (x-xc)*cos(alpha) - (y-yc)*sin(alpha) for x, y in zip(xcoords, ycoords)]
    yrot = [yc + (x-xc)*sin(alpha) + (y-yc)*cos(alpha) for x, y in zip(xcoords, ycoords)]

    return xrot, yrot


def CircleCurve(x0, R, theta, N):
    thetarad = pi * theta / 180.0
    pts = [Point(x0[0] + R*cos(a), x0[1] + R*sin(a), x0[2])
           for a in np.linspace(thetarad, thetarad + 2*pi, N+1)]

    return ip.CubicCurve(pts=pts, boundary=ip.PERIODIC)


def GenerateCircularDisc(x0, R, a, N):
    alpha = a * pi / 180

    thetas = np.linspace(0, 2*pi, N+1)
    xs = [x0[0] + R * cos(t + alpha) for t in thetas]
    ys = [x0[1] + R * sin(t + alpha) for t in thetas]
    zs = [x0[2] for _ in thetas]
    return ip.CubicCurve(x=xs, y=ys, z=zs, boundary=ip.PERIODIC).ReParametrize(0.0, 2*pi*R)


def GradedSpace(start, step, factor, N):
    def GeneratorGradedSpace(start, step, factor, N):
        for _ in xrange(N):
            yield start
            start += step
            step *= factor
    return list(GeneratorGradedSpace(float(start), float(step), float(factor), N))


def LinSpace(*args, **kwargs):
    return list(linspace(*args, **kwargs))


def MkInternalMesh(wv):
    pass


default_params = {
    'wingfile': 'NREL_5_MW.xinp',            # XML file for wing definition
    'out': 'NREL_wing_mesh_3D',              # Output filename
    'debug': False,                          # Debug mode produces intermediate g2-files
    'nprocs': 800,                           # Optimize for a given number of CPUs
    'order': 2,                              # 2 = linear, 3 = quadratic, etc. (what a dumb convention)

    'adds': 4,                               # Number of intermediate refinements
    'addpt': 4,                              # Airfoil index around which to perform intermediate refinement

    'gradedLen': True,                       # Whether to use graded refinement lengthwise
    'extrudeLen': 10.0,
    'dzTip': 0.006,                          # Element size in z-direction near tip
    'dzJoin': 0.2,                           # Element size in z-direction near join
    'NlenBase': 20,                          # Number of elements lengthwise at the base
    'Nlen': 70,                              # Number of elements for half a wing

    'bndlayer': 4,                           # Number of elements in boundary layer

    'te': 2e-2,                              # Size of trailing edge
    'teC_fac': 50,                           # Size of trailing edge for cylinders, in terms of `te`
    'R': 10.0,                               # Radius of the o-mesh, the sidelength of the full mesh
                                             # is four times as large

    'Nte': 9,                                # Cells for trailing edge
    'Nback': 27,                             # Cells for back part of wing
    'Nfront': 38,                            # Cells for front part of wing

    'Nrad': 70,                              # Number of elements radially in the o-mesh
    'NradSq': 8,                             # Number of elements radially outside the o-mesh
    'NradP': 4,                              # Number of patches radially

    'NlenP': 18,                             # Number of patches lengthwise

    'tip': True,                             # Whether to include tip

    'back': 4.0,

    'OpenFOAM': False,
}


if __name__ == '__main__':

    ParseArgs(sys.argv[1:], default_params)

    def fix(d):
        for k in ['z', 'theta', 'chord', 'ac', 'ao']:
            d[k] = float(d[k])
        return d

    wingdef_tree = ET.parse(default_params['wingfile'])
    default_params['wingdef'] = [namedtuple('Section', s.attrib.keys())(**fix(s.attrib))
                         for s in wingdef_tree.getroot()]
    default_params['teC'] = default_params['teC_fac'] * default_params['te']
    default_params['computed'] = {}

    params = namedtuple('Params', default_params.keys())(**default_params)
    PrepareOutput(params)

    # create non-dimensional sections
    ndsections = []                            # the normalized wing sections
    sections = []
    ste = []                                   # length of trailing edge modification
    sUpper = []                                # curve length for upper sections
    sTotal = []                                # total curve length for sections

    for secdef in params.wingdef:

        # Special treatment for cylinders
        if secdef.foil == 'cylinder':
            x0 = Point(secdef.chord * (0.25 - secdef.ao + secdef.ac), 0, secdef.z)
            s = GenerateCircularDisc(x0, 0.5*secdef.chord, secdef.theta, 200)
            knots = s.GetKnots()

            ndsections.append(s)
            ste.append(params.teC)
            sUpper.append(0.5 * knots[-1])
            sTotal.append(knots[-1])

            continue

        # Modification for trailing edge thickness
        x2D_nd, y2D_nd = TrailingEdgeModificationFt(secdef.foil, params.te / secdef.chord)

        # Modification for chord length
        x2Dc = [secdef.chord * (x - (secdef.ao + 0.25 - secdef.ac)) for x in x2D_nd]
        y2Dc = [secdef.chord * y for y in y2D_nd]

        # Rotate curve
        x0 = Point(0.0, 0.0, 0.0)
        x2D, y2D = RotateCoordinates(x2Dc, y2Dc, x0, secdef.theta)

        # Find tangent vector at trailing edge
        tau1 = Point(x2D[1]  - x2D[0],  y2D[1]  - y2D[0],  0); tau1.Normalize()
        tau2 = Point(x2D[-1] - x2D[-2], y2D[-1] - y2D[-2], 0); tau2.Normalize()

        # Generate points on wing section curve
        wingsection = ip.CubicCurve(x=x2D, y=y2D, z=[secdef.z]*len(x2D),
                                    boundary=ip.TANGENT, der=[tau1, tau2])
        knots = wingsection.GetKnots()

        WriteG2('out/wingsection.g2', [wingsection])

        # Trailing edge curve
        ptste, wingte = TrailingEdgeCurve(wingsection)
        imid = (len(ptste) - 1) / 2

        # Resample wing curve
        s1 = knots[(len(knots) - 1) / 2]
        s2 = knots[-1] - s1
        ds = abs(ptste[1] - ptste[0])
        r1 = ComputeFactor(0.5*s1, ds, 100)
        r2 = ComputeFactor(0.5*s2, ds, 100)

        ss =  GradedSpace(0,      ds,                r1, 101)
        ss += GradedSpace(ss[-1], ss[-1]-ss[-2], 1.0/r1, 101)[1:-1] + [s1]
        ss += GradedSpace(ss[-1], ds,                r2, 101)[1:]
        ss += GradedSpace(ss[-1], ss[-1]-ss[-2], 1.0/r2, 101)[1:-1] + [knots[-1]]

        pts1 = [wingsection.Evaluate(s) for s in ss]
        pts2 = ptste[imid:-1] + pts1 + ptste[1:imid+1]

        knots2 = CurveLengthParametrization(pts2, False)
        knotste = wingte.GetKnots()
        sUpper.append(knots2[(len(knots2) - 1) / 2])
        sTotal.append(knots2[-1])
        ste.append(knotste[-1])

        # Interpolate curve using Hermite interpolation
        ndsections.append(ip.CubicCurve(pts=pts2, boundary=ip.PERIODIC))

    if params.debug:
        WriteG2('out/ndsections.g2', ndsections)
        for i, c in enumerate(ndsections):
            WriteG2('out/ndsections/{:03}.g2'.format(i),
                    [c, NonRationalCurve(LineSegment(c.Evaluate(c.GetKnots()[0]),
                                                     c.EvaluateTangent(c.GetKnots()[0]),
                                                     relative=True))])


    # Loop over sections and sample mesh points
    for i, (sT, sU, st, nds) in enumerate(zip(sTotal, sUpper, ste, ndsections)):
        factor = 0.5
        s1 = 0.5 * st                     # Length of trailing edge modification
        ds1 = s1 / params.Nte                    # Grid spacing at trailing edge

        # If shit is fucked up near the hub, try changing this number first
        dx2 = 2.5e-3 if i >= 4 else 5e-2
        r2, r3, _ = GradingTwoSided(sU-s1, ds1, dx2, params.Nback-1, params.Nfront-1,
                                    1.1, 0.9, 1.0e-12, 200)
        r4, r5, _ = GradingTwoSided(sT-sU-s1, dx2, ds1, params.Nfront-1, params.Nback-1,
                                    1.1, 0.9, 1.0e-12, 200)

        knots  = list(np.linspace(0, s1, params.Nte+1))
        knots += GradedSpace(knots[-1], ds1, 1.0/r2, params.Nback+1)[1:]
        knots += GradedSpace(knots[-1], knots[-1]-knots[-2], r3, params.Nfront+1)[1:]
        knots += GradedSpace(knots[-1], knots[-1]-knots[-2], 1.0/r4, params.Nfront+1)[1:]
        knots += GradedSpace(knots[-1], knots[-1]-knots[-2], r5, params.Nback+1)[1:]
        knots += list(np.linspace(knots[-1], nds.GetKnots()[-1], params.Nte+1))[1:]

        pts = [nds.Evaluate(k) for k in knots]

        # Generate uniform parametrization
        sections.append(ip.CubicCurve(pts=pts, t=list(np.linspace(0, 1, len(pts))),
                                      boundary=ip.PERIODIC))


    # Intermediate interpolation to prevent crossing meshlines
    if params.adds > 0:

        def add_between(secs, i):
            pas, pbs = GetCurvePoints(secs[i-1]), GetCurvePoints(secs[i])
            pts = [(pa + pb)/2 for pa, pb in zip(pas, pbs)]
            return secs[:i] + [ip.CubicCurve(pts=pts, t=list(np.linspace(0, 1, len(pts))),
                                             boundary=ip.PERIODIC)] + secs[i:]

        for i in xrange(params.adds):
            sections = add_between(sections, params.addpt+i)
        for _ in xrange(params.adds):
            sections = add_between(sections, params.addpt+params.adds+1)
        sections = sections[:params.addpt+params.adds] + sections[params.addpt+params.adds+1:]

    if params.debug:
        WriteG2('out/sections.g2', sections)
        for i, c in enumerate(sections):
            WriteG2('out/sections/{:03}.g2'.format(i),
                    [c, NonRationalCurve(LineSegment(c.Evaluate(c.GetKnots()[0]),
                                                     c.EvaluateTangent(c.GetKnots()[0]).Normalize() * 0.1,
                                                     relative=True))])

    if len(params.wingdef) > 1:
        wing1 = LoftCurves(sections, order=4)

        if params.debug:
            WriteG2('out/wing1.g2', wing1)

        knots1, knots2 = wing1.GetKnots(True)

        zcoord = [wing1.Evaluate(knots1[0], k2)[2] for k2 in knots2]
        zcoord[1:4] = repeat(zcoord[0], 3)
        zcoord[-4:-1] = repeat(zcoord[-1], 3)

        coeffs = list(wing1)

        # Curve length parametrization
        p1, p2 = wing1.GetOrder()
        wing = Surface(p1, p2, knots1, zcoord, coeffs, False)

        if params.debug:
            WriteG2('out/wing.g2', wing)

        wingsecs = ResampleWing(wing, params)

    else:
        wingsecs = sections


    if params.debug:
        WriteG2('out/wingsecs.g2', wingsecs)
        for i, c in enumerate(wingsecs):
            WriteG2('out/wingsecs/{:03}.g2'.format(i),
                    [c, NonRationalCurve(LineSegment(c.Evaluate(c.GetKnots()[0]),
                                                     c.EvaluateTangent(c.GetKnots()[0]),
                                                     relative=True))])


    # Generate circular sections
    if len(params.wingdef) > 1:
        thetapts = [Point(0, 0, th) for th in map(attrgetter('theta'), params.wingdef)]
        thetacurve = ip.CubicCurve(pts=thetapts, t=map(attrgetter('z'), params.wingdef), boundary=ip.NATURAL)
        theta = lambda z: thetacurve.Evaluate(z)[2]
    else:
        theta = lambda _: params.wingdef[0].theta

    Ntot = 2 * (params.Nte + params.Nback + params.Nfront)
    circlesecs = []
    for ws in wingsecs:
        pt = ws.Evaluate(ws.GetKnots()[0])
        x0 = Point(0, 0, pt[2])
        angle = theta(pt[2])
        circlesecs.append(CircleCurve(x0, params.R, angle, Ntot))

    if params.debug:
        WriteG2('out/circlesecs.g2', circlesecs)
        for i, c in enumerate(circlesecs):
            WriteG2('out/circlesecs/{:03}.g2'.format(i), c)





    # Generate 2D cross sections
    crossecs = []
    normalf = lambda _, crv, kts: [crv.EvaluateTangent(k).Rotate(Point(0,0,1), pi/2).Normalize()
                                   for k in kts]
    init_normalf = lambda crv, kts: normalf(None, crv, kts)
    for i, (ws, cs) in enumerate(zip(wingsecs, circlesecs)):
        sys.stdout.write('\rTFI: %i/%i...' % (i+1, len(wingsecs))); sys.stdout.flush()

        ws.FlipParametrization()
        knots = ws.GetKnots()
        x1  = ws.Evaluate(knots[0])
        t1  = ws.EvaluateTangent(knots[0]).Normalize()
        n1 = Point(-t1[1], t1[0], 0)

        cs.FlipParametrization()
        knots = cs.GetKnots()
        x2 = cs.Evaluate(knots[0])
        t2 = cs.EvaluateTangent(knots[0]).Normalize()
        n2 = Point(-t2[1], t2[0], 0)

        c31 = ip.CubicCurve(pts=[x1, x2], boundary=ip.TANGENT, der=[n1, n2])

        s3 = c31.GetKnots()[-1]
        r3 = 1.0/0.9
        coeff3 = (1.0 - pow(r3,params.Nrad)) / (1.0-r3)
        knots3 = GradedSpace(0.0, s3/coeff3, r3, params.Nrad) + [s3]

        points3 = [c31.Evaluate(k) for k in knots3]
        c3 = ip.CubicCurve(pts=points3, t=knots3, boundary=ip.TANGENT, der=[n1, n2])
        c3.InsertKnot(0.5 * (knots3[0] + knots3[1]))

        if params.debug:
            WriteG2('out/tfi/{:03}.g2'.format(i),
                    [ws, cs, c3,
                     NonRationalCurve(LineSegment(x1, t1, relative=True)),
                     NonRationalCurve(LineSegment(x1, n1, relative=True))])

        crossecs.append(tfi.OrthogonalSurface([ws, cs, c3, c3], init_normalf, normalf,
                                              ranges=[(0,15), (-15,0)], fac_blend=0.9,
                                              fac_smooth=0.98, nsweeps=1, bnd_layer=params.bndlayer))

    sys.stdout.write('\n'); sys.stdout.flush()

    if params.debug:
        WriteG2('out/crossecs.g2', crossecs)
        for i, c in enumerate(crossecs):
            WriteG2('out/crossecs/{:03}.g2'.format(i), c)





    # Split each section into 8
    print 'Splitting patches in angular direction...'

    ku, _ = crossecs[-1].GetKnots()
    nq = (len(ku) - 1) / 4
    ns = nq / 2
    nl = nq - ns
    idx_split = [0, ns, nq, nq+nl, 2*nq, 2*nq+ns, 3*nq, 3*nq+nl, 4*nq]

    crossecs8, square8, innersecs, outersecs = [], [], [], []
    r = 1.0 / 0.9
    for csec in crossecs:
        knotsu, knotsv = csec.GetKnots()
        nu = len(knotsu)

        pt = csec.Evaluate(knotsu[0], knotsv[0])
        z = pt[2]

        splitKnots = [knotsu[i] for i in idx_split]
        splitPts = [csec.Evaluate(k, knotsv[-1]) for k in splitKnots]

        nv1, nv2 = knotsv[0], knotsv[-1]
        pt1 = csec.Evaluate(knotsu[0], knotsv[-2])
        pt2 = csec.Evaluate(knotsu[0], knotsv[-1])

        # Points on outer square
        points = [Point( 2.0*params.R,  0.0,           z),
                  Point( 2.0*params.R, -2.0*params.R,  z),
                  Point( 0.0,          -2.0*params.R,  z),
                  Point(-2.0*params.R, -2.0*params.R,  z),
                  Point(-2.0*params.R,  0.0,           z),
                  Point(-2.0*params.R,  2.0*params.R,  z),
                  Point( 0.0,           2.0*params.R,  z),
                  Point( 2.0*params.R,  2.0*params.R,  z),
                  Point( 2.0*params.R,  0.0,           z)]

        slen = min([abs(splitpt - pt) for splitpt, pt in zip(splitPts, points)[:-1]])
        dy = abs(pt1 - pt2)

        for sKnp, sKnn, sPtp, sPtn, ptp, ptn in zip(splitKnots[:-1],  splitKnots[1:],
                                                    splitPts[:-1],    splitPts[1:],
                                                    points[:-1],      points[1:]):
            p1, p2 = Point(sKnp, nv1), Point(sKnn, nv2)
            innersec = csec.GetSubSurf(p1, p2)

            edges = innersec.GetEdges()
            c1 = edges[2]

            knots2 = edges[2].GetKnots()
            pts = [(1.0-s)*ptp + s*ptn for s in np.linspace(0, 1, len(knots2))]
            c2 = ip.LinearCurve(pts=pts).RaiseOrder(2)

            c3 = ip.InterpolateCurve(np.matrix([sPtp, ptp]), [0, 1], [0, 0, 1, 1])
            c4 = ip.InterpolateCurve(np.matrix([sPtn, ptn]), [0, 1], [0, 0, 1, 1])
            rf = ComputeFactor(slen, dy, params.NradSq)
            sst = GradedSpace(0.0, dy/slen, rf, params.NradSq)[1:]
            for c in [c3, c4]:
                c.RaiseOrder(2)
                for k in sst:
                    c.InsertKnot(k)

            outersec = tfi.LinearSurface([c1, c2, c3, c4])

            innersecs.append(innersec)
            outersecs.append(outersec)
            crossecs8.append(innersec)
            crossecs8.append(outersec)

    grouper = lambda n, it, fv: list(izip_longest(*[iter(it)]*n, fillvalue=fv))

    if params.debug:
        WriteG2('out/inner.g2', innersecs)
        for i, c in enumerate(grouper(8, innersecs, None)):
            WriteG2('out/inner/{:03}.g2'.format(i), c)
        WriteG2('out/outer.g2', outersecs)
        for i, c in enumerate(grouper(8, outersecs, None)):
            WriteG2('out/outer/{:03}.g2'.format(i), c)





    if params.tip:
        print 'Generating tip...'

        wv = []
        for secs in [innersecs[i::8] for i in range(8)][::-1]:
            wv.append(LoftSurfaces(secs[-5:], range(5), params.order))

        cs = [o.Clone() for o in crossecs8[-16::2]]
        os = [o.Clone() for o in outersecs[-8:]]

        tip = WingTip(wv, params)
        # sphere = OMeshTip(params.wingdef[-1].z, params.R, cs, tip)
        # tfi_surfaces = TipTFISurfaces(cs, tip, sphere, params)
        # tipvols = TipCircleVolumes(cs, tip, sphere, tfi_surfaces)
        # outervols = OuterVolumes(os, tipvols, params)

        # if params.debug:
        #     WriteG2('out/tipvols.g2', tipvols + outervols)




    if len(params.wingdef) == 1:
        srfs = []

        # Kill me now
        for i, o in zip(innersecs, outersecs):
            ikus, ikvs = i.GetKnots()
            okus, okvs = o.GetKnots()

            eikus = [ikus[0], (ikus[0]+ikus[1])/2] + ikus[1:-1] + [(ikus[-2]+ikus[-1])/2, ikus[-1]]
            eokus = [okus[0], (okus[0]+okus[1])/2] + okus[1:-1] + [(okus[-2]+okus[-1])/2, okus[-1]]

            eikvs = [ikvs[0], (ikvs[0]+ikvs[1])/2] + ikvs[1:]
            eokvs = okvs[1:-1] + [(okvs[-2]+okvs[-1])/2, okvs[-1]]

            pts = []
            for ikv in eikvs:
                for iku in eikus:
                    pts.append(i.Evaluate(iku, ikv))
            for okv in eokvs:
                for oku in eokus:
                    pts.append(o.Evaluate(oku, okv))

            ukts = map(float, range(len(ikus)))
            vkts = map(float, range(len(ikvs) + len(okvs) - 1))
            eukts = [ukts[0], (ukts[0]+ukts[1])/2] + ukts[1:-1] + [(ukts[-2]+ukts[-1])/2, ukts[-1]]
            evkts = [vkts[0], (vkts[0]+vkts[1])/2] + vkts[1:-1] + [(vkts[-2]+vkts[-1])/2, vkts[-1]]
            ukts = [ukts[0]]*3 + ukts + [ukts[-1]]*3
            vkts = [vkts[0]]*3 + vkts + [vkts[-1]]*3

            srf = ip.InterpolateSurface(pts, eukts, evkts, ukts, vkts)

            srfs.append(srf)

        vols_master = [ExtrudeSurface(s, Point(0, 0, 1), params.extrudeLen) for s in srfs]
        for v in vols_master:
            v.RaiseOrder(0, 0, 2)
            UniformVolume(v, 3, params.NlenBase + 2*params.Nlen)
            v.SwapParametrization(1, 2)
            v.FlipParametrization(0)

    else:

        print 'Constructing master volumes...'

        innervols_master = []
        outervols_master = []
        vols_master = []

        if not params.tip:
            tipvols = [None]*12
            outervols = [None]*12

        for secs, tv in zip([innersecs[i::8] for i in range(8)], tipvols[:8]):
            if tv:
                tvsecs = [tv.GetConstParSurf(k, 1).LowerOrder(2,2) for k in tv.GetKnots()[1][1:]]
            else:
                tvsecs = []
            allsecs = [s.ReParametrize() for s in secs + tvsecs]
            innervols_master.append(LoftSurfaces(allsecs, range(len(allsecs)), order=4))

        for secs, tv in zip([outersecs[i::8] for i in range(8)], outervols[:8]):
            if tv:
                tvsecs = [tv.GetConstParSurf(k, 0).SwapParametrization() for k in tv.GetKnots()[0][1:]]
            else:
                tvsecs = []
            allsecs = [s.ReParametrize() for s in secs + tvsecs]
            outervols_master.append(LoftSurfaces(allsecs, range(len(allsecs)), order=4))

        for iv, ov in zip(innervols_master, outervols_master):
            iss = [iv.GetConstParSurf(k, 1) for k in iv.GetKnots()[1]]
            iss += [ov.GetConstParSurf(k, 1) for k in ov.GetKnots()[1][1:]]
            vols_master.append(LoftSurfaces(iss, range(len(iss)), order=4))

        if params.tip:
            for iv, ov, flip in zip(tipvols[8:], outervols[8:], [False, True, False, True]):
                ivs = [iv.GetConstParSurf(k, 2).LowerOrder(2,2) for k in iv.GetKnots()[2]]
                ovs = [ov.GetConstParSurf(k, 0).SwapParametrization().FlipParametrization(1)
                       for k in ov.GetKnots()[0][1:]]
                allsecs = ivs + ovs
                vols_master.append(LoftSurfaces(allsecs, range(len(allsecs)), order=4))


        if params.debug:
            WriteG2('out/innervols_master.g2', innervols_master)
            WriteG2('out/outervols_master.g2', outervols_master)

    if params.debug:
        WriteG2('out/full_master.g2', vols_master)





    print 'Subdividing...'

    radials = []
    for v in vols_master:
        ku, kv, kw = v.GetKnots()

        ret = []
        N = len(kw) - 1
        for i in range(params.NradP):
            k1 = kw[i * N / params.NradP]
            k2 = kw[(i+1) * N / params.NradP]
            ret.append(v.GetSubVol([ku[0], kv[0], k1], [ku[-1], kv[-1], k2]))
        radials.append(ret)

    if params.debug:
        WriteG2('out/radials.g2', list(chain.from_iterable(radials)))

    length_vols = []
    for sector in radials[:8]:
        final_sector = []
        for column in sector:
            final_column = []
            ku, kv, kw = column.GetKnots()

            N = len(kv) - 1
            for i in range(params.NlenP):
                k1 = kv[i * N / params.NlenP]
                k2 = kv[(i+1) * N / params.NlenP]
                final_column.append(column.GetSubVol([ku[0], k1, kw[0]], [ku[-1], k2, kw[-1]]))
            final_sector.append(final_column)
        length_vols.append(final_sector)

    out_vols = list(chain.from_iterable(chain.from_iterable(length_vols)))
    if params.tip:
        out_vols += list(chain.from_iterable(radials[:8]))





    if params.back > 0:
        backvols = []
        for l in length_vols[0][-1] + length_vols[7][-1]:
            srf = l.GetFaces()[5]
            backvols.append(ExtrudeSurface(srf, Point(1, 0, 0), params.R * params.back))
        for v in backvols:
            v.RaiseOrder(0, 0, 2)
            UniformVolume(v, 3, int(ceil(params.back * params.NradSq)))
        out_vols += backvols

    if params.debug:
        WriteG2('out/vols.g2', out_vols)





    maxdist, npts = 0.0, 0
    N = 50
    for ray_a, ray_b in zip(length_vols, length_vols[1:] + [length_vols[0]]):
        for col_a, col_b in zip(ray_a, ray_b):
            for va, vb in zip(col_a, col_b):
                a_ukts, a_vkts, a_wkts = va.GetKnots(with_multiplicities=True)
                b_ukts, b_vkts, b_wkts = vb.GetKnots(with_multiplicities=True)
                assert(a_wkts == b_wkts)
                assert(a_vkts == b_vkts)

                a_ukts, a_vkts, a_wkts = va.GetKnots()
                b_ukts, b_vkts, b_wkts = vb.GetKnots()
                assert(a_wkts == b_wkts)
                assert(a_vkts == b_vkts)

                vtest = np.linspace(a_vkts[0], a_vkts[-1], N)
                wtest = np.linspace(a_wkts[0], a_wkts[-1], N)

                a_pts = [va.Evaluate(a_ukts[0], v, w) for v in vtest for w in wtest]
                b_pts = [vb.Evaluate(b_ukts[-1], v, w) for v in vtest for w in wtest]

                npts += len(a_pts)
                maxdist = max(max(abs(a-b) for a, b in zip(a_pts, b_pts)), maxdist)

    for ray in length_vols:
        for col_i, col_o in zip(ray[:-1], ray[1:]):
            for va, vb in zip(col_i, col_o):
                a_ukts, a_vkts, a_wkts = va.GetKnots(with_multiplicities=True)
                b_ukts, b_vkts, b_wkts = vb.GetKnots(with_multiplicities=True)
                assert(a_ukts == b_ukts)
                assert(a_vkts == b_vkts)

                a_ukts, a_vkts, a_wkts = va.GetKnots()
                b_ukts, b_vkts, b_wkts = vb.GetKnots()
                assert(a_ukts == b_ukts)
                assert(a_vkts == b_vkts)

                utest = np.linspace(a_ukts[0], a_ukts[-1], N)
                vtest = np.linspace(a_vkts[0], a_vkts[-1], N)

                a_pts = [va.Evaluate(u, v, a_wkts[-1]) for u in utest for v in vtest]
                b_pts = [vb.Evaluate(u, v, b_wkts[0]) for u in utest for v in vtest]

                npts += len(a_pts)
                maxdist = max(max(abs(a-b) for a, b in zip(a_pts, b_pts)), maxdist)

    for ray in length_vols:
        for col in ray:
            for va, vb in zip(col[:-1], col[1:]):
                a_ukts, a_vkts, a_wkts = va.GetKnots(with_multiplicities=True)
                b_ukts, b_vkts, b_wkts = vb.GetKnots(with_multiplicities=True)
                assert(a_ukts == b_ukts)
                assert(a_wkts == b_wkts)

                a_ukts, a_vkts, a_wkts = va.GetKnots()
                b_ukts, b_vkts, b_wkts = vb.GetKnots()
                assert(a_ukts == b_ukts)
                assert(a_wkts == b_wkts)

                utest = np.linspace(a_ukts[0], a_ukts[-1], N)
                wtest = np.linspace(a_wkts[0], a_wkts[-1], N)

                a_pts = [va.Evaluate(u, a_vkts[-1], w) for u in utest for w in wtest]
                b_pts = [vb.Evaluate(u, b_vkts[0], w) for u in utest for w in wtest]

                npts += len(a_pts)
                maxdist = max(max(abs(a-b) for a, b in zip(a_pts, b_pts)), maxdist)

    print 'Maximal node gap found: %s (tested %i points)' % (maxdist, npts)





    for v in out_vols:
        v.LowerOrder(4 - params.order, 4 - params.order, 4 - params.order)

    numberer = Numberer()
    numberer.AddPatches(out_vols)

    numberer.AddGroup('inflow', 'volume', length_vols[3][-1] + length_vols[4][-1])
    numberer.AddGroup('out_left', 'volume', length_vols[0][-1])
    numberer.AddGroup('out_right', 'volume', length_vols[7][-1])
    numberer.AddGroup('slip_left', 'volume', length_vols[1][-1] + length_vols[2][-1])
    numberer.AddGroup('slip_right', 'volume', length_vols[5][-1] + length_vols[6][-1])
    numberer.AddGroup('slip_left_in', 'volume', length_vols[2][-1])
    numberer.AddGroup('slip_right_in', 'volume', length_vols[5][-1])
    numberer.AddGroup('inner', 'volume', list(chain.from_iterable([length_vols[k][0] for k in xrange(8)])))
    numberer.AddGroup('btm', 'volume', [v[0] for v in chain.from_iterable(length_vols)])
    numberer.AddGroup('top', 'volume', [v[-1] for v in chain.from_iterable(length_vols)])

    if params.back > 0:
        numberer.AddGroup('outflow', 'volume', backvols)
        numberer.AddGroup('btm', 'volume', [backvols[0], backvols[len(backvols)/2]])
        numberer.AddGroup('top', 'volume', [backvols[len(backvols)/2-1], backvols[-1]])
        numberer.AddGroup('slip_left_back', 'volume', backvols[:len(backvols)/2])
        numberer.AddGroup('slip_right_back', 'volume', backvols[len(backvols)/2:])
    else:
        numberer.AddGroup('outflow', 'volume', length_vols[0][-1] + length_vols[7][-1])

    numberer.AddBoundary('hub', [('btm', 'face', [2])])
    numberer.AddBoundary('antihub', [('top', 'face', [3])])
    numberer.AddBoundary('inflow', [('inflow', 'face', [5]),
                                    ('slip_left_in', 'edge', [10]),
                                    ('slip_right_in', 'edge', [11])])
    numberer.AddBoundary('outflow', [('outflow', 'face', [5])])
    numberer.AddBoundary('slipwall_left', [('slip_left', 'face', [5]),
                                           ('out_left', 'edge', [10])])
    numberer.AddBoundary('slipwall_right', [('slip_right', 'face', [5]),
                                            ('out_right', 'edge', [11])])
    numberer.AddBoundary('wing', [('inner', 'face', [4])])

    if params.back > 0:
        numberer.AddBoundary('slipwall_left', [('slip_left_back', 'face', [0]),
                                               ('slip_left_back', 'edge', [10])])
        numberer.AddBoundary('slipwall_right', [('slip_right_back', 'face', [1]),
                                                ('slip_right_back', 'edge', [11])])


    numberer.Renumber(params.nprocs)


    if params.debug:
        for g in numberer.Groups():
            numberer.WriteGroup(g, 'out/group-%s.g2' % g)
        for g in numberer.Boundaries():
            numberer.WriteBoundary(g, 'out/boundary-%s.g2' % g)


    numberer.WriteEverything(params.out)
    numberer.PrintLoadBalance()

    if params.OpenFOAM:
        print 'Converting to OpenFOAM format...'
        f = InputFile('%s.xinp' % params.out)
        f.writeOpenFOAM(params.out)
