from utils import *
from slice import Slice

class VolumetricSlice(Slice):
    """This class represents a three-dimensional slice of the blade geometry. It is subclassed from
    Slice since much of the functionality is shared."""

    def __init__(self, params):
        """Private constructor. Foreign code should use `from_slices` or `from_slice`."""
        self.p = params

    @classmethod
    def from_slices(cls, slices):
        """Creates a volumetric slice by lofting slices."""
        new = cls(slices[0].p)
        for attr in slices[0].components():
            setattr(new, attr, deep_loft([getattr(s, attr) for s in slices]))
        return new

    @classmethod
    def from_slice(cls, s):
        """Creates a volumetric slice by extruding a single slice by unit distance."""
        new = cls(s.p)
        for attr in s.components():
            setattr(new, attr, deep_extrude(getattr(s, attr)))
        return new

    @classmethod
    def from_attrs(cls, params, **kwargs):
        """Crates a volumetric slice with the given components. Foreign code should not call this
        constructor."""
        new = cls(params)
        for attr, val in kwargs.iteritems():
            setattr(new, attr, val)
        return new

    def z(self):
        """Prevents foreign code from relying on Slice.z()."""
        raise NotImplementedError("VolumetricSlice.z() doesn't make sense")

    def subdivide(self):
        """Subdivides a volumetric slice lengthwise. This function shadows Slice.subdivide(). One
        should call Slice.subdivide() before constructing a volumetric slice."""
        slices = [{} for _ in xrange(self.p.p_length)]

        for attr in self.components():
            temp = deep_subdivide(getattr(self, attr), self.p.p_length, 2)
            for i, s in enumerate(slices):
                s[attr] = deep_index(temp, i)

        return [VolumetricSlice.from_attrs(self.p, **s) for s in slices]

    def push_topsets(self, n, complete=False):
        """Adds the boundaries in this slice to the numberer (wing, inflow, outflow and
        slipwalls). Also adds hub and antihub if complete is True."""
        self._bnd_wing(n, edge_kind='face')

        # The Slice._bnd_xyz methods will work for volumes if we provide the right input
        kwargs = {'edge_kind': 'face', 'vx_kind': 'edge', 'vx_add': 8}
        self._bnd_slipwall(n, **kwargs)
        self._bnd_flow(n, 'inflow', **kwargs)
        self._bnd_flow(n, 'outflow', **kwargs)

        if complete:
            self.push_hub(n)
            self.push_antihub(n)

        self._grp_rigid(n, kind='volume', edge_kind='face', vx_kind='edge', vx_add=8)

    def push_hub(self, n):
        """Adds the hub boundary to the numberer."""
        self._bnd_hub(n, 'hub')

    def push_antihub(self, n):
        """Adds the antihub boundary to the numberer."""
        self._bnd_hub(n, 'antihub')

    def _bnd_hub(self, n, hub):
        """Adds hub or antihub boundary to the numberer."""
        # Index of face (up or down)
        face = 4 if hub == 'hub' else 5

        # The easy part
        for attr in self.components():
            n.AddTopologySet(hub, (flatten_objects(getattr(self, attr)), 'face', face))

        # Utility function for adding edges
        edge_add = 2 if hub == 'antihub' else 0
        def edge(patches, target, idx):
            if patches:
                n.AddTopologySet(target, (patches, 'edge', idx + edge_add))

        # Add edges on the wing interface to the correct boundary
        target = hub if getattr(self.p, 'wing_' + hub) == hub else 'wing'
        patches = [q[0] for q in self.inner_left + self.inner_right]
        n.AddTopologySet(target, (patches, 'edge', edge_add))

        # Add edges on the inflow interface to the correct boundary
        target = hub if getattr(self.p, 'in_' + hub) == hub else 'inflow'
        patches = []
        if self.p.ext_a: patches += [q[-1] for q in self.ahead]
        if self.p.ext_as: patches += [q[-1] for c in self.corners_ahead for q in c]
        if self.p.ext_not_a: patches += [self.inner_left[3][-1], self.inner_right[3][-1]]
        if self.p.ext_s_not_a: patches += self.left[-1] + self.right[-1]
        edge(patches, target, 1)

        # Add edges on the outflow interface to the correct boundary
        target = hub if getattr(self.p, 'out_' + hub) == hub else 'outflow'
        patches = []
        if self.p.ext_b: patches += [q[0] for q in self.behind]
        if self.p.ext_bs: patches += [q[0] for q in c for c in self.corners_behind]
        if self.p.ext_not_b: edge([self.inner_left[0][-1], self.inner_right[0][-1]], target, 1)
        if self.p.ext_s_not_b: patches += self.left[0] + self.right[0]
        edge(patches, target, 0)

        # Add edges on the slipwall interfaces to the correct boundaries
        target_l = hub if getattr(self.p, 'slip_' + hub) == hub else 'slipwall_left'
        target_r = hub if getattr(self.p, 'slip_' + hub) == hub else 'slipwall_right'
        patches_l, patches_r = [], []

        if self.p.ext_s:
            patches_l += [q[0] for q in self.left]
            patches_r += [q[-1] for q in self.right]
        else:
            edge([q[-1] for q in self.inner_left[1:3]], target_l, 1)
            edge([q[-1] for q in self.inner_right[1:3]], target_r, 1)

        for attr in ['behind', 'ahead']:
            if getattr(self.p, 'ext_s' + attr[0]):
                patches_l += getattr(self, 'corners_' + attr)[0][0]
                patches_r += getattr(self, 'corners_' + attr)[1][-1]
            elif getattr(self.p, 'ext_' + attr[0] + '_not_s'):
                patches_l += getattr(self, attr)[0]
                patches_r += getattr(self, attr)[1]

        edge(patches_l, target_l, 4)
        edge(patches_r, target_r, 5)
