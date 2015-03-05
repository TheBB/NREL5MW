import numpy as np
import os

from math import pi, sin, cos, sqrt

from GoTools import Point, Curve, Surface, WriteG2
from GoTools.CurveFactory import Circle, IntersectCurve, LineSegment, NonRationalCurve
from GoTools.VolumeFactory import ExtrudeSurface
from GeoUtils.CurveUtils import CurveLengthParametrization, GetCurvePoints
from GeoUtils.Factory import LoftBetween
from GeoUtils.Refinement import UniformCurve, GeometricRefineCurve, GeometricRefineSurface
import GeoUtils.Interpolate as ip
from GeoUtils.IO import Numberer, InputFile
import GeoUtils.TFI as tfi

from utils import *


COMPONENTS = ['inner_left', 'inner_right', 'behind', 'ahead', 'left', 'right',
              'corners_behind', 'corners_ahead']


class AirFoil(object):

    def __init__(self, filled=False, volumetric=False, params=None):
        self.filled = filled
        self.volumetric = volumetric
        if params:
            self.p = params


    @classmethod
    def from_wd(cls, params, wd):
        obj = cls(filled=False, volumetric=False)
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
        obj = cls(filled=False, volumetric=False)
        obj._from_pts(pts)
        obj.theta = theta

        return obj


    @classmethod
    def average(cls, afa, afb):
        assert(afa.curve.GetKnots() == afb.curve.GetKnots())

        pts = [(pta + ptb) / 2 for pta, ptb in zip(GetCurvePoints(afa.curve),
                                                   GetCurvePoints(afb.curve))]
        return cls.from_pts(pts, (afa.theta + afb.theta) / 2)


    @classmethod
    def loft_volumetric(cls, airfoils):
        new = cls(filled=True, volumetric=True, params=airfoils[0].p)

        for attr in COMPONENTS:
            if hasattr(airfoils[0], attr):
                setattr(new, attr, deep_loft([getattr(af, attr) for af in airfoils]))

        return new


    def dummy_volumetric(self):
        def deep_extrude(objs):
            if type(objs) is Surface:
                vol = ExtrudeSurface(objs, ez, 1.0)
                return vol
            else:
                return [deep_extrude(q) for q in objs]

        for attr in COMPONENTS:
            if hasattr(self, attr):
                setattr(self, attr, deep_extrude(getattr(self, attr)))

        self.volumetric = True


    def translate(self, pt):
        assert(self.filled)

        new = AirFoil(filled=True, volumetric=False, params=self.p)
        new.curve = deep_translate(self.curve, pt)

        for attr in COMPONENTS:
            if hasattr(self, attr):
                setattr(new, attr, deep_translate(getattr(self, attr), pt))

        return new


    def subdivide_volumetric(self, n):
        new = [AirFoil(filled=True, volumetric=True, params=self.p)
               for _ in xrange(n)]

        for attr in COMPONENTS:
            if hasattr(self, attr):
                temp = deep_subdivide(getattr(self, attr), n, 2)
                for i, af in enumerate(new):
                    setattr(af, attr, deep_index(temp, i))

        return new


    def objects(self):
        if not self.filled and not self.volumetric:
            return [self.curve]

        objects = []
        for attr in COMPONENTS:
            if hasattr(self, attr):
                objects.extend(flatten_objects(getattr(self, attr)))

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

        for p in split:
            p.FlipParametrization(0)

        self.inner_left = split[4:][::-1]
        self.inner_right = split[:4]
        self.filled = True


    def fill_sides(self):
        ext = lambda ps, direc, dist, n: [extend(p.GetEdges()[2] if type(p) is Surface else p,
                                                 direc, dist, n) for p in ps]

        if self.p.sides > 0:
            self.left  = ext(self.inner_left[1:3], -ey, self.p.sides, self.p.n_sides)
            self.right = ext(self.inner_right[1:3], ey, self.p.sides, self.p.n_sides)
            deep_orient(self.left, 'flipv', 'swap')
            deep_orient(self.right, 'flipu', 'swap')

        if self.p.behind > 0:
            edges = [self.inner_left[0], self.inner_right[0]]
            self.behind = ext(edges, ex, self.p.behind, self.p.n_behind)
            deep_orient(self.behind, 'flipu', 'flipv')
            if self.p.sides > 0:
                edges = [self.left[0].GetEdges()[0], self.right[0].GetEdges()[0]]
                self.corners_behind = ext(edges, ex, self.p.behind, self.p.n_behind)
                deep_orient(self.corners_behind, 'flipv')

        if self.p.ahead > 0:
            edges = [self.inner_left[-1], self.inner_right[-1]]
            self.ahead = ext(edges, -ex, self.p.ahead, self.p.n_ahead)
            if self.p.sides > 0:
                edges = [self.left[-1].GetEdges()[2], self.right[-1].GetEdges()[2]]
                self.corners_ahead = ext(edges, -ex, self.p.ahead, self.p.n_ahead)


    def subdivide(self):
        self.inner_left = deep_subdivide(self.inner_left, self.p.p_inner, 1)
        self.inner_right = deep_subdivide(self.inner_right, self.p.p_inner, 1)

        if self.p.sides > 0:
            self.left = deep_subdivide(self.left, self.p.p_sides, 0)
            self.right = deep_subdivide(self.right, self.p.p_sides, 0)

        if self.p.behind > 0:
            self.behind = deep_subdivide(self.behind, self.p.p_behind, 1)
            if self.p.sides > 0:
                temp = deep_subdivide(self.corners_behiond, self.p.p_sides, 0)
                self.corners_behind = deep_subdivide(temp, self.p.p_behind, 1)

        if self.p.ahead > 0:
            self.ahead = subdivide(self.ahead, self.p.p_ahead, 1)
            if self.p.sides > 0:
                temp = deep_subdivide(self.corners_ahead, self.p.p_sides, 0)
                self.corners_ahead = deep_subdivide(temp, self.p.p_ahead, 1)


    def lower_order(self, target, objs=None):
        def lower(objs):
            if type(objs) in [Surface, Volume]:
                orders = [order - target for order in objs.GetOrder()]
                objs.LowerOrder(*orders)
            elif type(objs) is list:
                for obj in objs:
                    lower(obj)

        for attr in COMPONENTS:
            if hasattr(self, attr):
                lower(getattr(self, attr))


    def output(self, path):
        n = Numberer()
        self.put_to_numberer(n)

        if self.volumetric:
            self._bnd_hub(n, 'hub')
            self._bnd_hub(n, 'antihub')

        if self.p.debug:
            n.WriteBoundaries(path)

        if self.p.walldistance:
            n.AddWallGroup('wing')

        # Final output
        def progress(s):
            sys.stdout.write('\r' + s)
            sys.stdout.flush()

        n.Renumber(self.p.nprocs)
        n.WriteEverything(path, display=False)

        if self.p.format == 'OpenFOAM':
            f = InputFile('%s.xinp' % path)
            f.writeOpenFOAM(os.path.dirname(path))

            for postfix in ['.xinp', '.g2', '_nodenumbers.hdf5', '_nodenumbers.xml']:
                os.remove(path + postfix)


    def put_to_numberer(self, n):
        for attr in COMPONENTS:
            if hasattr(self, attr):
                n.AddPatches(flatten_objects(getattr(self, attr)))

        self._bnd_wing(n)
        self._bnd_slipwall(n)
        self._bnd_flow(n, 'inflow')
        self._bnd_flow(n, 'outflow')


    def hub_to_numberer(self, n):
        self._bnd_hub(n, 'hub')


    def antihub_to_numberer(self, n):
        self._bnd_hub(n, 'antihub')


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


    def _bnd_wing(self, n):
        patches = [q[0] for q in self.inner_left + self.inner_right]
        kind = 'face' if self.volumetric else 'edge'
        n.AddBoundary('wing', (patches, kind, 2))


    def _bnd_slipwall(self, n):
        edge_kind = 'face' if self.volumetric else 'edge'
        vx_kind = 'edge' if self.volumetric else 'vertex'
        vx_add = 8 if self.volumetric else 0

        def edge(side, patches, idx):
            n.AddBoundary('slipwall_%s' % side, (patches, edge_kind, idx))

        def vx(side, patches, idx):
            n.AddBoundary('slipwall_%s' % side, (patches, vx_kind, idx + vx_add))

        # Add the easy part
        if self.p.sides > 0:
            edge('left', [q[0] for q in self.left], 0)
            edge('right', [q[-1] for q in self.right], 1)
        else:
            edge('left', [q[-1] for q in self.inner_left[1:3]], 3)
            edge('right', [q[-1] for q in self.inner_right[1:3]], 3)

        # Add the remainder of the edges, possibly corner patches
        for e in ['behind', 'ahead']:
            if getattr(self.p, e) == 0:
                continue
            if self.p.sides > 0:
                corners = getattr(self, 'corners_' + e)
                l, r = corners[0][0], corners[-1][-1]
            else:
                patches = getattr(self, e)
                l, r = patches[0], patches[-1]
            edge('left', l, 0)
            edge('right', r, 1)

        # In case of behind/ahead and no sides, there are some vertices in the middle of the slipwall
        if self.p.ahead > 0 and self.p.sides == 0:
            vx('left', [self.inner_left[-1][-1]], 2)
            vx('right', [self.inner_right[-1][-1]], 3)
        if self.p.behind > 0 and self.p.sides == 0:
            vx('left', [self.inner_left[0][-1]], 3)
            vx('right', [self.inner_right[0][-1]], 2)

        # In case of open inflow, ensure the corner vertices belong to the slipwall
        if self.p.in_slip == 'slip':
            if self.p.ahead == 0 and self.p.sides == 0:
                vx('left', self.inner_left[3][-1], 2)
                vx('right', self.inner_right[3][-1], 3)
            else:
                if self.p.ahead > 0 and self.p.sides == 0:
                    l, r = self.ahead[0][-1], self.ahead[1][-1]
                elif self.p.ahead == 0 and self.p.sides > 0:
                    l, r = self.left[1][0], self.right[1][-1]
                elif self.p.ahead > 0 and self.p.sides > 0:
                    l, r = self.corners_ahead[0][0][-1], self.corners_ahead[1][-1][-1]
                vx('left', l, 2)
                vx('right', r, 3)

        # In case of open outflow, ensure the corner vertices belong to the slipwall
        if self.p.out_slip == 'slip':
            if self.p.behind == 0 and self.p.sides == 0:
                vx('left', self.inner_left[0][-1], 3)
                vx('right', self.inner_right[0][-1], 2)
            else:
                if self.p.behind > 0 and self.p.sides == 0:
                    l, r = self.behind[0][0], self.behind[1][0]
                elif self.p.behind == 0 and self.p.sides > 0:
                    l, r = self.left[0][0], self.right[0][-1]
                elif self.p.behind > 0 and self.p.sides > 0:
                    l, r = self.corners_behind[0][0][0], self.corners_behind[1][-1][0]
                vx('left', l, 0)
                vx('right', r, 1)


    def _bnd_flow(self, n, name):
        out = name == 'outflow'

        ext_name = 'behind' if out else 'ahead'      # Name of extension property
        v_idx = 0 if out else -1                     # Index of patches in the v-direction
        v_edge = 2 if out else 3                     # Patch-index of edge on appropriately oriented patch
        l_vx = 0 if out else 2                       # Patch-index of left vertex on A.O.P.
        i_fidx = 0 if out else 3                     # Index of far patch in inner arrays
        i_nidx = 1 if out else 2                     # Index of near patch in inner arrays
        i_or = 2 if out else 3                       # Patch-index of the vertex closest to the middle
                                                     # on the left side of an inner patch

        edge_kind = 'face' if self.volumetric else 'edge'
        vx_kind = 'edge' if self.volumetric else 'vertex'
        vx_add = 8 if self.volumetric else 0

        def edge(patches, idx):
            n.AddBoundary(name, (patches, edge_kind, idx))

        def vx(patches, idx):
            n.AddBoundary(name, (patches, vx_kind, idx + vx_add))
    
        closed = getattr(self.p, name[:-4] + '_slip') != 'slip'
        extended = getattr(self.p, ext_name)
        if extended:
            ext_patches = getattr(self, ext_name)
        if extended and self.p.sides > 0:
            corners = getattr(self, 'corners_' + ext_name)

        j_or = 5 - i_or

        # Add the easy part
        if extended:
            edge([q[v_idx] for q in ext_patches], v_edge)
        else:
            edge([self.inner_left[i_fidx][-1], self.inner_right[i_fidx][-1]], 3)

        # Add the remainder of the edges, possibly corner patches
        if self.p.sides > 0:
            if extended:
                p = [q[v_idx] for c in corners for q in c]
            else:
                p = self.left[v_idx] + self.right[v_idx]
            edge(p, v_edge)

        # In case of sides and not extended, there are some vertices in the middle of the flow
        if self.p.sides > 0 and not extended:
            vx(self.inner_left[i_nidx][-1], i_or)
            vx(self.inner_right[i_nidx][-1], j_or)

        # In case of closed flow, ensure the corner vertices belong to the flow
        if closed:
            if not extended and self.p.sides == 0:
                vx(self.inner_left[i_nidx][-1], i_or)
                vx(self.inner_right[i_nidx][-1], j_or)
            else:
                if extended and self.p.sides == 0:
                    l, r = ext_patches[0][v_idx], ext_patches[1][v_idx]
                elif not extended and self.p.sides > 0:
                    l, r = self.left[v_idx][0], self.right[v_idx][-1]
                elif extended and self.p.sides > 0:
                    l, r = corners[0][0][v_idx], corners[1][-1][v_idx]
                vx(l, l_vx)
                vx(r, l_vx + 1)


    def _bnd_hub(self, n, hub):
        face = 4 if hub == 'hub' else 5

        # The easy part
        for attr in COMPONENTS:
            if hasattr(self, attr):
                n.AddBoundary(hub, (flatten_objects(getattr(self, attr)), 'face', face))

        edge_add = 2 if hub == 'antihub' else 0
        def edge(patches, target, idx):
            if patches:
                n.AddBoundary(target, (patches, 'edge', idx + edge_add))

        # Add edges on the inflow interface to the correct boundary
        target = hub if getattr(self.p, 'in_' + hub) == hub else 'inflow'
        if self.p.ahead > 0:
            patches = [q[-1] for q in self.ahead]
            if self.p.sides > 0:
                patches += [q[-1] for q in c for c in self.corners_ahead]
        else:
            patches = [self.inner_left[3][-1], self.inner_right[3][-1]]
            if self.p.sides > 0:
                patches += self.left[-1] + self.right[-1]
        edge(patches, target, 1)

        # Add edges on the outflow interface to the correct boundary
        target = hub if getattr(self.p, 'out_' + hub) == hub else 'outflow'
        patches = []
        if self.p.behind > 0:
            patches += [q[0] for q in self.behind]
            if self.p.sides > 0:
                patches += [q[0] for q in c for c in self.corners_behind]
        else:
            edge([self.inner_left[0][-1], self.inner_right[0][-1]], target, 1)
            if self.p.sides > 0:
                patches += self.left[0] + self.right[0]
        edge(patches, target, 0)

        # # Add edges on the slipwall interfaces to the correct boundaries
        target_l = hub if getattr(self.p, 'slip_' + hub) == hub else 'slipwall_left'
        target_r = hub if getattr(self.p, 'slip_' + hub) == hub else 'slipwall_right'
        patches_l, patches_r = [], []
        if self.p.sides > 0:
            patches_l += [q[0] for q in self.left]
            patches_r += [q[-1] for q in self.right]
            for attr in ['behind', 'ahead']:
                if getattr(self.p, attr) > 0:
                    patches_l += getattr(self, 'corners_' + attr)[0][0]
                    patches_r += getattr(self, 'corners_' + attr)[1][-1]
        else:
            edge([q[-1] for q in self.inner_left[1:3]], target_l, 1)
            edge([q[-1] for q in self.inner_right[1:3]], target_r, 1)
            for attr in ['behind', 'ahead']:
                if getattr(self.p, attr) > 0:
                    patches_l += getattr(self, attr)[0]
                    patches_r += getattr(self, attr)[1]
        edge(patches_l, target_l, 4)
        edge(patches_r, target_r, 5)
