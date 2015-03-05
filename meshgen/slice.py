from GoTools.CurveFactory import Circle, IntersectCurve, LineSegment
from GeoUtils.IO import Numberer, InputFile
from GeoUtils.Refinement import GeometricRefineCurve, GeometricRefineSurface, UniformCurve
import GeoUtils.TFI as tfi

from airfoil import AirFoil
from utils import *


COMPONENTS = {'inner_left', 'inner_right',
              'behind', 'ahead', 'left', 'right',
              'corners_behind', 'corners_ahead'}


class Slice(object):

    def __init__(self, af):
        self.p = af.p

        trailing = self._make_trailing(af.curve)
        inner = self._make_inner(trailing, af.curve)
        split_inner = self._split_inner(inner)
        split_outer = self._make_outer(split_inner)
        split = [merge_surfaces(i, m, 1) for i, m in zip(split_inner, split_outer)]

        deep_orient(split, 'flipu')

        self.inner_left = split[4:][::-1]
        self.inner_right = split[:4]

    def components(self):
        return set(self.__dict__.keys()) & COMPONENTS

    def objects(self):
        return chain.from_iterable(flatten_objects(getattr(self, attr))
                                   for attr in self.components())

    def z(self):
        return self.inner_left[0][0][0][2]

    def extend(self):
        ext = lambda ps, direc, dist, n: [extend(p.GetEdges()[2] if type(p) is Surface else p,
                                                 direc, dist, n) for p in ps]

        if self.p.ext_s:
            self.left  = ext(self.inner_left[1:3], -ey, self.p.sides, self.p.n_sides)
            self.right = ext(self.inner_right[1:3], ey, self.p.sides, self.p.n_sides)
            deep_orient(self.left, 'flipv', 'swap')
            deep_orient(self.right, 'flipu', 'swap')

        if self.p.ext_b:
            edges = [self.inner_left[0], self.inner_right[0]]
            self.behind = ext(edges, ex, self.p.behind, self.p.n_behind)
            deep_orient(self.behind, 'flipu', 'flipv')

        if self.p.ext_bs:
            edges = [self.left[0].GetEdges()[0], self.right[0].GetEdges()[0]]
            self.corners_behind = ext(edges, ex, self.p.behind, self.p.n_behind)
            deep_orient(self.corners_behind, 'flipv')

        if self.p.ext_a:
            edges = [self.inner_left[-1], self.inner_right[-1]]
            self.ahead = ext(edges, -ex, self.p.ahead, self.p.n_ahead)

        if self.p.ext_as:
            edges = [self.left[-1].GetEdges()[2], self.right[-1].GetEdges()[2]]
            self.corners_ahead = ext(edges, -ex, self.p.ahead, self.p.n_ahead)

    def subdivide(self):
        SPEC = [(['inner_left', 'inner_right'], ('p_inner', 1)),
                (['left', 'right'], ('p_sides', 0)),
                ('behind', ('p_behind', 1)),
                ('ahead', ('p_ahead', 1)),
                ('corners_behind', [('p_sides', 0), ('p_behind', 1)]),
                ('corners_ahead', [('p_sides', 0), ('p_ahead', 1)])]

        for attrs, divs in SPEC:
            if type(attrs) is str:
                attrs = [attrs]
            if type(divs) is tuple:
                divs = [divs]

            for attr in set(attrs) & self.components():
                temp = getattr(self, attr)
                for num, direction in divs:
                    temp = deep_subdivide(temp, getattr(self.p, num), direction)
                setattr(self, attr, temp)

    def lower_order(self):
        for attr in self.components():
            deep_lower_order(getattr(self, attr), self.p.order)

    def output(self, path):
        n = Numberer()
        self.push_patches(n)
        self.push_boundaries(n, complete=True)

        if self.p.debug:
            n.WriteBoundaries(path)

        if self.p.walldistance:
            n.AddWallGroup('wing')

        n.Renumber(self.p.nprocs)
        n.WriteEverything(path, display=False)

        if self.p.format == 'OpenFOAM':
            convert_openfoam(path)

    def push_patches(self, n):
        for attr in self.components():
            n.AddPatches(flatten_objects(getattr(self, attr)))

    def push_boundaries(self, n, complete=False):
        self._bnd_wing(n)
        self._bnd_slipwall(n)
        self._bnd_flow(n, 'inflow')
        self._bnd_flow(n, 'outflow')

    def _bnd_wing(self, n, kind='edge', idx=2):
        patches = [q[0] for q in self.inner_left + self.inner_right]
        n.AddBoundary('wing', (patches, kind, idx))

    def _bnd_slipwall(self, n, edge_kind='edge', vx_kind='vertex', vx_add=0):
        def edge(side, patches, idx):
            n.AddBoundary('slipwall_%s' % side, (patches, edge_kind, idx))
        def vx(side, patches, idx):
            n.AddBoundary('slipwall_%s' % side, (patches, vx_kind, idx + vx_add))

        # Add the easy part
        if self.p.ext_s:
            edge('left', [q[0] for q in self.left], 0)
            edge('right', [q[-1] for q in self.right], 1)
        else:
            edge('left', [q[-1] for q in self.inner_left[1:3]], 3)
            edge('right', [q[-1] for q in self.inner_right[1:3]], 3)

        # Add the remainder of the edges, possibly corner patches
        for e in {'behind', 'ahead'} & self.components():
            if self.p.ext_s:
                corners = getattr(self, 'corners_' + e)
                l, r = corners[0][0], corners[-1][-1]
            else:
                patches = getattr(self, e)
                l, r = patches[0], patches[-1]
            edge('left', l, 0)
            edge('right', r, 1)

        # In case of behind/ahead and no sides, there are some vertices in the middle of the slipwall
        if self.p.ext_a and not self.p.ext_s:
            vx('left', [self.inner_left[-1][-1]], 2)
            vx('right', [self.inner_right[-1][-1]], 3)
        if self.p.ext_b and not self.p.ext_s:
            vx('left', [self.inner_left[0][-1]], 3)
            vx('right', [self.inner_right[0][-1]], 2)

        # In case of open inflow, ensure the corner vertices belong to the slipwall
        if self.p.in_slip == 'slip':
            if self.p.ext_not_as:
                l, r = self.inner_left[3][-1], self.inner_right[3][-1]
            elif self.p.ext_a_not_s:
                l, r = self.ahead[0][-1], self.ahead[1][-1]
            elif self.p.ext_s_not_a:
                l, r = self.left[1][0], self.right[1][-1]
            elif self.p.ext_as:
                l, r = self.corners_ahead[0][0][-1], self.corners_ahead[1][-1][-1]
            vx('left', l, 2)
            vx('right', r, 3)
                
        # In case of open outflow, ensure the corner vertices belong to the slipwall
        if self.p.out_slip == 'slip':
            if self.p.ext_bs:
                vx('left', self.inner_left[0][-1], 3)
                vx('right', self.inner_right[0][-1], 2)
            else:
                if self.p.ext_b_not_s:
                    l, r = self.behind[0][0], self.behind[1][0]
                elif self.p.ext_s_not_b:
                    l, r = self.left[0][0], self.right[0][-1]
                elif self.p.ext_bs:
                    l, r = self.corners_behind[0][0][0], self.corners_behind[1][-1][0]
                vx('left', l, 0)
                vx('right', r, 1)

    def _bnd_flow(self, n, flow, edge_kind='edge', vx_kind='vertex', vx_add=0):
        out = flow == 'outflow'

        ext_name = 'behind' if out else 'ahead'      # Name of extension property
        v_idx = 0 if out else -1                     # Index of patches in the v-direction
        v_edge = 2 if out else 3                     # Patch-index of edge on appropriately oriented patch
        l_vx = 0 if out else 2                       # Patch-index of left vertex on A.O.P.
        i_fidx = 0 if out else 3                     # Index of far patch in inner arrays
        i_nidx = 1 if out else 2                     # Index of near patch in inner arrays
        i_or = 2 if out else 3                       # Patch-index of the vertex closest to the middle
                                                     # on the left side of an inner patch
        j_or = 5 - i_or

        def edge(patches, idx):
            n.AddBoundary(flow, (patches, edge_kind, idx))
        def vx(patches, idx):
            n.AddBoundary(flow, (patches, vx_kind, idx + vx_add))

        closed = getattr(self.p, flow[:-4] + '_slip') != 'slip'
        extended = getattr(self.p, 'ext_' + ext_name[0])
        if extended:
            ext_patches = getattr(self, ext_name)
        if extended and self.p.ext_s:
            corners = getattr(self, 'corners_' + ext_name)

        # Add the easy part
        if extended:
            edge([q[v_idx] for q in ext_patches], v_edge)
        else:
            edge([self.inner_left[i_fidx][-1], self.inner_right[i_fidx][-1]], 3)

        # Add the remainder of the edges, possibly corner patches
        if self.p.ext_s:
            if extended:
                edge([q[v_idx] for c in corners for q in c], v_edge)
            else:
                edge(self.left[v_idx] + self.right[v_idx], v_edge)

        # In case of sides and not extended, there are some vertices in the middle of the flow
        if self.p.ext_s and not extended:
            vx(self.inner_left[i_nidx][-1], i_or)
            vx(self.inner_right[i_nidx][-1], j_or)
        
        # In case of closed flow, ensure the corner vertices belong to the flow
        if closed:
            if not extended and not self.p.ext_s:
                vx(self.inner_left[i_nidx][-1], i_or)
                vx(self.inner_right[i_nidx][-1], j_or)
            else:
                ext_es = (extended, self.p.ext_s)
                l, r = {(True,  False): (ext_patches[0][v_idx], ext_patches[1][v_idx]),
                        (False, True):  (self.left[v_idx][0], self.right[v_idx][-1]),
                        (True,  True):  (corners[0][0][v_idx], corners[1][-1][v_idx])}[ext_es]
                vx(l, l_vx)
                vx(r, l_vx + 1)

    def _make_trailing(self, curve):
        z = curve[0][2]
        k = curve.GetKnots()[0]
        p_inner = curve.Evaluate(k)
        n_inner = curve.EvaluateTangent(k).Rotate(ez, -pi/2).Normalize()

        p_outer = Point(2 * self.p.radius, 0, z)

        trailing = ip.CubicCurve(pts=[p_inner, p_outer], der=[n_inner, ex],
                                 boundary=ip.TANGENT)

        circle = Circle(ez * z, self.p.radius, ez)
        point = IntersectCurve(trailing, circle)[1][0]
        knot = trailing.GetParameterAtPoint(point)[0]

        trailing.InsertKnot(knot)
        trailing = trailing.GetSubCurve(0, knot)
        fac = self.p.radial_grading(knot)
        GeometricRefineCurve(trailing, fac, self.p.n_circle - 1)

        self.inner_fac = fac

        return trailing

    def _make_inner(self, trailing, curve):
        z = trailing[0][2]
        point = trailing.Evaluate(trailing.GetKnots()[-1])
        theta = np.arctan(point[1] / point[0]) * 180 / pi
        radius = abs(trailing.Evaluate(trailing.GetKnots()[-1]) - ez * z)

        circle =  mkcircle(ez * z, radius, theta,
                           len(curve.GetKnots()) - 1)

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
            [curve, circle, trailing, trailing],
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

    def _make_outer(self, inners):
        z = inners[0][0][2]
        radius = 2 * self.p.radius

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
