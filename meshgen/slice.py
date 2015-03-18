from itertools import *

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
    """This class represents a two-dimensional slice of the blade geometry."""

    def __init__(self, params):
        """Private constructor. Foreign code should call `from_airfoil`."""
        self.p = params

    @classmethod
    def from_airfoil(cls, airfoil):
        """Constructs a slice from an airfoil. The parameters are inherited from the given
        airfoil. The constructor will make all the TFI computations necessary to produce the
        surrounding O-mesh."""
        new = cls(airfoil.p)

        trailing, fac = new._make_trailing(airfoil.curve)
        inner = new._make_inner(trailing, airfoil.curve)
        split_inner = new._split_inner(inner)
        split_outer = new._make_outer(split_inner, fac)
        split = [merge_surfaces(i, m, 1) for i, m in zip(split_inner, split_outer)]

        # Ensure right-handedness
        deep_orient(split, 'flipu')

        new.inner_left = split[4:][::-1]
        new.inner_right = split[:4]

        return new

    def translate(self, pt):
        """Create a new slice by translating this slice along a given vector."""
        new = Slice(self.p)
        for attr in self.components():
            setattr(new, attr, deep_translate(getattr(self, attr), pt))
        return new

    def components(self):
        """Return the names of all the patch components in this slice, a (possibly proper) subset of
        `COMPONENTS`."""
        return set(self.__dict__.keys()) & COMPONENTS

    def objects(self):
        """Required for debug output."""
        return chain.from_iterable(flatten_objects(getattr(self, attr))
                                   for attr in self.components())

    def z(self):
        """Returns the z-location of this slice."""
        return self.inner_left[0][0][0][2]

    def extend(self):
        """Creates extension patches according to the given parameters (behind, ahead and side)."""

        # Utility function for extending a list of surfaces or curves
        ext = lambda ps, direc, dist, n: [extend(p.GetEdges()[2] if type(p) is Surface else p,
                                                 direc, dist, n) for p in ps]

        # Sides
        if self.p.ext_s:
            self.left  = ext(self.inner_left[1:3], -ey, self.p.sides, self.p.n_sides)
            self.right = ext(self.inner_right[1:3], ey, self.p.sides, self.p.n_sides)
            deep_orient(self.left, 'flipv', 'swap')
            deep_orient(self.right, 'flipu', 'swap')

        # Behind
        if self.p.ext_b:
            edges = [self.inner_left[0], self.inner_right[0]]
            self.behind = ext(edges, ex, self.p.behind, self.p.n_behind)
            deep_orient(self.behind, 'flipu', 'flipv')

        # Corners behind
        if self.p.ext_bs:
            edges = [self.left[0].GetEdges()[0], self.right[0].GetEdges()[0]]
            self.corners_behind = ext(edges, ex, self.p.behind, self.p.n_behind)
            deep_orient(self.corners_behind, 'flipv')

        # Ahead
        if self.p.ext_a:
            edges = [self.inner_left[-1], self.inner_right[-1]]
            self.ahead = ext(edges, -ex, self.p.ahead, self.p.n_ahead)

        # Corners ahead
        if self.p.ext_as:
            edges = [self.left[-1].GetEdges()[2], self.right[-1].GetEdges()[2]]
            self.corners_ahead = ext(edges, -ex, self.p.ahead, self.p.n_ahead)

    def subdivide(self):
        """Subdivides the components according to the given parameters."""
        # A list of all subdivisions to perform
        # Each entry is tuple with:
        # - The name of list of names of component(s) to subdivide
        # - The subdivision or list of subdivisions to perform, in order.
        # Each subdivision is a tuple with:
        # - The name of a parameter giving the number of patches
        # - The direction in which to do the subdivision
        SPEC = [(['inner_left', 'inner_right'], ('p_inner', 1)),
                (['left', 'right'], ('p_sides', 0)),
                ('behind', ('p_behind', 1)),
                ('ahead', ('p_ahead', 1)),
                ('corners_behind', [('p_sides', 0), ('p_behind', 1)]),
                ('corners_ahead', [('p_sides', 0), ('p_ahead', 1)])]

        # Actually perform the subdivisions
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
        """Lowers the order of all the components according to the parameters."""
        for attr in self.components():
            deep_lower_order(getattr(self, attr), self.p.order)

    def output(self, path, custom=None):
        """Outputs this slice by itself to the given path."""
        n = Numberer()
        self.push_patches(n)
        self.push_topsets(n, complete=True)

        if self.p.debug:
            n.WriteTopologySets(path)

        if self.p.walldistance:
            n.AddWallGroup('wing')

        n.Renumber(self.p.nprocs)
        n.WriteEverything(path, display=False)

        if self.p.format == 'OpenFOAM':
            convert_openfoam(path)
        if self.p.format == 'IFEM':
            self.p.postprocess_xinp(custom)

    def push_patches(self, n):
        """Adds all the patches to the numberer."""
        for attr in self.components():
            n.AddPatches(flatten_objects(getattr(self, attr)))

    def push_topsets(self, n, complete=False):
        """Adds the topology sets in this slice to the numberer (wing, inflow, outflow, slipwalls
        and rigid). Call `push_patches` first."""
        self._bnd_wing(n)
        self._bnd_slipwall(n)
        self._bnd_flow(n, 'inflow')
        self._bnd_flow(n, 'outflow')
        self._grp_rigid(n)

    def _grp_rigid(self, n, kind='face', edge_kind='edge', vx_kind='vertex', vx_add=0):
        if self.p.p_rigid == 0:
            return

        def edge(patches, idx):
            n.AddTopologySet('rigid', (patches, edge_kind, idx))
        def vx(patches, idx):
            n.AddTopologySet('rigid', (patches, vx_kind, idx + vx_add))

        # Volumetric part
        patches = [q[:self.p.p_rigid] for q in self.inner_left + self.inner_right]
        patches = list(chain.from_iterable(patches))
        n.AddTopologySet('rigid', (patches, kind, []))

        # Boundary part
        if self.p.p_rigid < self.p.p_inner:
            edge([q[self.p.p_rigid] for q in self.inner_left + self.inner_right], 2)
        else:
            if self.p.ext_s:
                edge([q[-1] for q in self.left], 1)
                edge([q[0] for q in self.right], 0)
            if self.p.ext_b:
                edge([q[-1] for q in self.behind], 3)
            if self.p.ext_a:
                edge([q[0] for q in self.ahead], 2)
            if self.p.ext_bs:
                vx(self.corners_behind[0][-1][-1], 3)
                vx(self.corners_behind[1][0][-1], 2)
            if self.p.ext_as:
                vx(self.corners_ahead[0][-1][0], 1)
                vx(self.corners_ahead[1][0][0], 0)

    def _bnd_wing(self, n, edge_kind='edge'):
        """Adds the wing boundary to the numberer."""
        patches = [q[0] for q in self.inner_left + self.inner_right]
        n.AddTopologySet('wing', (patches, edge_kind, 2))

    def _bnd_slipwall(self, n, edge_kind='edge', vx_kind='vertex', vx_add=0):
        """Adds the slipwall boundaries to the numberer."""
        # Utility functions for adding an edge or a vertex
        def edge(side, patches, idx):
            n.AddTopologySet('slipwall_%s' % side, (patches, edge_kind, idx))
        def vx(side, patches, idx):
            n.AddTopologySet('slipwall_%s' % side, (patches, vx_kind, idx + vx_add))

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
        if self.p.ext_a_not_s:
            vx('left', [self.inner_left[-1][-1]], 2)
            vx('right', [self.inner_right[-1][-1]], 3)
        if self.p.ext_b_not_s:
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
        """Adds the in- or outflow boundary to the numberer."""
        # True if outflow, false if inflow
        out = flow == 'outflow'

        # All the things that are different between the outflow and the inflow
        ext_name = 'behind' if out else 'ahead'      # Name of extension property
        v_idx = 0 if out else -1                     # Index of patches in the v-direction
        v_edge = 2 if out else 3                     # Patch-index of edge on appropriately oriented patch
        l_vx = 0 if out else 2                       # Patch-index of left vertex on A.O.P.
        i_fidx = 0 if out else 3                     # Index of far patch in inner arrays
        i_nidx = 1 if out else 2                     # Index of near patch in inner arrays
        i_or = 2 if out else 3                       # Patch-index of the vertex closest to the middle
                                                     # on the left side of an inner patch
        j_or = 5 - i_or

        # Utility functions for adding edges and vertices
        def edge(patches, idx):
            n.AddTopologySet(flow, (patches, edge_kind, idx))
        def vx(patches, idx):
            n.AddTopologySet(flow, (patches, vx_kind, idx + vx_add))

        # Whether this boundary is closed, or has an extension
        closed = getattr(self.p, flow[:-4] + '_slip') != 'slip'
        extended = getattr(self.p, 'ext_' + ext_name[0])

        # Grab the extension patches and the corner patches, if needed
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
        """Returns a curve from the trailing edge to the O-mesh circle. Also returns the geometric
        grading factor used to refine it."""
        # Get the z-value, the point and the normal at the trailing edge
        z = curve[0][2]
        k = curve.GetKnots()[0]
        p_inner = curve.Evaluate(k)
        n_inner = curve.EvaluateTangent(k).Rotate(ez, -pi/2).Normalize()

        # The corresponding point on the square (not the circle)
        p_outer = Point(2 * self.p.radius, 0, z)

        # A suitable curve connecting the two
        trailing = ip.CubicCurve(pts=[p_inner, p_outer], der=[n_inner, ex],
                                 boundary=ip.TANGENT)

        # Find the intersection point on the circle
        circle = Circle(ez * z, self.p.radius, ez)
        point = IntersectCurve(trailing, circle)[1][0]
        knot = trailing.GetParameterAtPoint(point)[0]

        # Get the subcurve and geometrically refine
        trailing.InsertKnot(knot)
        trailing = trailing.GetSubCurve(0, knot)
        fac = self.p.radial_grading(knot)
        GeometricRefineCurve(trailing, fac, self.p.n_circle - 1)

        return trailing, fac

    def _make_inner(self, trailing, curve):
        """Creates the inner part of the O-mesh (the circle)."""
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
        """Splits the inner O-mesh in eight patches."""
        ksu, ksv = inner.GetKnots()
        ksu_split = [ksu[i] for i in self.p.angular_splits()]

        splits = []
        for kp, kn in zip(ksu_split[:-1], ksu_split[1:]):
            splits.append(inner.GetSubSurf((kp, ksv[0]), (kn, ksv[-1])))

        return splits

    def _make_outer(self, inners, inner_fac):
        """Creates the outer part of the O-mesh (the square)."""
        z = inners[0][0][2]
        radius = 2 * self.p.radius

        # Connecting points on the square
        outer_pts = [Point(radius, 0, z),
                     Point(radius, radius, z),
                     Point(0, radius, z),
                     Point(-radius, radius, z),
                     Point(-radius, 0, z),
                     Point(-radius, -radius, z),
                     Point(0, -radius, z),
                     Point(radius, -radius, z),
                     Point(radius, 0, z)]

        # Grading of the space inbetween
        kus, kvs = inners[0].GetKnots()
        ipt = inners[0].Evaluate(kus[0], kvs[-1])
        length = abs(ipt - outer_pts[0])
        ds = abs(ipt - inners[0].Evaluate(kus[0], kvs[-2])) * inner_fac
        fac = grading(length, ds, self.p.n_square)

        outers = []
        for opp, opn, inner in zip(outer_pts[:-1], outer_pts[1:], inners):
            kus, kvs = inner.GetKnots()

            # Create the outer curve
            temp_curve = LineSegment(opp, opn)
            UniformCurve(temp_curve, len(kus) - 2)
            out = ip.CubicCurve(pts=GetCurvePoints(temp_curve), boundary=ip.NATURAL, t=kus)

            # Perform the lofting and refinement
            surface = LoftBetween(inner.GetEdges()[2], out).RaiseOrder(0, 2)
            GeometricRefineSurface(surface, 2, fac, self.p.n_square - 1)

            # Useful for merging later
            diff = surface.GetKnots()[1][1]
            surface.ReParametrize(kus[0], kus[-1], kvs[-1], kvs[-1] + 1/diff)

            outers.append(surface)

        return outers
