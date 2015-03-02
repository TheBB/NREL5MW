# -*- coding: utf-8 -*-

import os
import sys
import xml.etree.ElementTree as xml

from itertools import product
from collections import namedtuple
from datetime import datetime
from math import sqrt
from shutil import rmtree
import subprocess
import yaml

from GoTools import WriteG2
from GeoUtils.IO import ParseArgs

from utils import grading, gradspace


defaults = {
    # Input and output
    'wingfile': 'wingdefs/NREL_5_MW.xinp',
    'out': 'NREL',
    'nprocs': 8,

    # Mesh generator
    'debug': False,
    'nprocs_mg': 4,
    'walldistance': False,
    'mesh_mode': '3d',

    # Global mesh parameters
    'order': 2,

    # Boundary conditions
    'in_slip': 'in',
    'in_hub': 'in',
    'in_antihub': 'in',
    'out_slip': 'slip',
    'out_hub': 'hub',
    'out_antihub': 'antihub',
    'slip_hub': 'hub',
    'slip_antihub': 'antihub',

    # Trailing edge
    'len_te': 2e-2,
    'len_te_cyl_fac': 50.0,

    # Angular resolution
    'n_te': 10,
    'n_back': 29,
    'n_front': 41,

    # Lengthwise resolution
    'length_mode': 'triple',
    'length': 0.0,
    'join_adds': 4,
    'join_index': 4,
    'n_base': 20,
    'n_length': 160,
    'd_join': 0.2,
    'd_tip': 0.006,

    # Radial resolution
    'radius': 10.0,
    'Re': 10000.0,
    'n_bndlayer': 4,
    'n_circle': 72,
    'n_square': 10,
    'smoothing': True,

    # Sides
    'behind': 4.0,
    'ahead': 0.0,
    'sides': 0.0,
    'n_behind': 36,
    'n_ahead': 2,
    'n_sides': 2,

    # Patches
    'p_inner': 2,
    'p_behind': 1,
    'p_ahead': 1,
    'p_sides': 1,
    'p_length': 20,
}


def usage():
    print """The following parameters are recognized:

    INPUT AND OUTPUT
    - wingfile: Path to XML file for wing definition
    - out: Path for writing output (will be deleted if exists)
    - nprocs: Number of processors to optimize for

    MESH GENERATOR
    - debug: Set to true for outputting intermediate geometry
    - nprocs_mg: Number of processors to use in mesh generation
    - walldistance: Set to true to compute wall distances
    - mesh_mode: Which kind of mesh to produce
      * 2d: 2D mesh of a single airfoil
        Make sure join_adds is zero, any other length parameter will be ignored
      * semi3d: Layered 2D meshes with a beam
      * 3d: 3D model with cut-off tip

    GLOBAL MESH PARAMETERS
    - order: Spline geometry order (2, 3 or 4)

    BOUNDARY CONDITIONS
    Each xxx_yyy determines whether the interface between boundaries xxx and yyy should
    belong to xxx or yyy.  Those are the only valid values.
    - in_slip: Inflow and slipwall
    - in_hub: Inflow and hub
    - in_antihub: Inflow and antihub
    - out_slip: Outflow and slipwall
    - out_hub: Inflow and hub
    - out_antihub: Inflow and antihub
    - slip_hub: Slipwall and hub
    - slip_antihub: Slipwall and antihub

    TRAILING EDGE
    - len_te: Size of the trailing edge modification
    - len_te_cyl_fac: How large should the “trailing edge” be in cylindrical
      airfoils (in terms of len_te)

    ANGULAR RESOLUTION
    - n_te: Number of elements for the trailing edge
    - n_back: Number of elements for the back of the airfoil
    - n_front: Number of elements for the front of the airfoil

    LENGTHWISE RESOLUTION
    - length_mode: How to distribute the airfoils in the z-direction
      * extruded: Only one airfoil should be specificed in the wing definition file
      * uniform: Uniformly distributed in the z-direction
      * double: Double-sided geometrically distributed airfoils, by giving the element
        size at the root and at the tip
      * triple: Triple-sided geometrically distributed airfoils, by giving the element
        size at the join (see below) and at the tip
    - length: Length of wing (in case of length mode extrude)
    - join_adds: Number of intermediate linear interpolation steps to perform at the
      join. The join is any airfoil with a sharp transition requiring this step to avoid
      self intersection. Set this to zero to disable
    - join_index: Index of the airfoil at the join, if any
    - n_base: Number of elements in the base (in case of length mode triple)
    - n_length: Number of lengthwise elements in total
    - d_join: Element size at the join or base (in case of length mode double or triple)
    - d_tip: Element size at the tip (in case of length mode double or triple)

    RADIAL RESOLUTION
    - radius: Radius of O-mesh
    - Re: Reynold's number
    - n_bndlayer: Number of elements in the boundary layer
    - n_circle: Number of elements in the O-mesh
    - n_square: Number of elements outside the O-mesh in the square
    - smoothing: True to turn on Laplacian smoothing during TFI

    EXTENSION PATCHES
    - behind: Distance to extend the outflow, in terms of radius
    - ahead: Distance to extend the inflow, in terms of radius
    - sides: Distance to extend the slipwalls, in terms of radius
    - n_behind: Number of elements behind
    - n_ahead: Number of elements ahead
    - n_sides: Number of elements on the sides

    SUBDIVISION
    - p_inner: Number of patches radially in the square
    - p_behind: Number of patches behind
    - p_ahead: Number of patches ahead
    - p_sides: Number of patches on the sides
    - p_length: Number of patches lengthwise"""


def fix_floats(dct, keys=['z', 'theta', 'chord', 'ac', 'ao']):
    for k in keys:
        dct[k] = float(dct[k])
    return dct


class Colors(object):
    WARNING = '\033[93m\033[1m'
    ERROR = '\033[91m\033[1m'
    END = '\033[0m'


def check_warn(test, msg):
    if not test:
        print Colors.WARNING + "WARNING: " + Colors.END + msg


def check_error(test, msg):
    if not test:
        print Colors.ERROR + "ERROR: " + Colors.END + msg
        sys.exit(0)


Section = namedtuple('Section', ['z', 'theta', 'chord', 'ac', 'ao', 'foil'])


class Params(object):

    def __init__(self, args=[]):
        ParseArgs(args, defaults)
        self.original = defaults
        self.__dict__.update(defaults)

        wingdef = xml.parse(self.wingfile)
        self.wingdef = [Section(**fix_floats(s.attrib)) for s in wingdef.getroot()]

        self.len_te_cyl = self.len_te_cyl_fac * self.len_te
        self.len_char = min(wd.chord for wd in self.wingdef)
        self.len_bndlayer = self.len_char / sqrt(self.Re)
        self.n_bndlayers = float(self.n_circle) / self.n_bndlayer

        self.behind *= self.radius
        self.ahead *= self.radius
        self.sides *= self.radius

        self.sanity_check()

        rmtree('out', ignore_errors=True)
        self.make_folder('out')

        rmtree(self.out, ignore_errors=True)
        self.make_folder(self.out)

        if self.mesh_mode == '2d':
            s = '2D mode'
        elif self.mesh_mode == 'semi3d':
            s = 'Semi3D mode (%i planes)' % (self.n_length + 1)
        elif self.mesh_mode == '3d':
            s = '3D mode'
        s += ' -- ' + {2: 'linear', 3: 'quadratic', 4: 'cubic'}[self.order] + ' geometry'
        print s


    def out_yaml(self):
        filename = os.path.join(self.out, 'parameters.yaml')

        time = datetime.now().isoformat()
        proc = subprocess.Popen(['git', '-C', os.path.abspath(os.path.dirname(__file__)),
                                 'rev-parse', 'HEAD'], stdout=subprocess.PIPE)
        commit, _ = proc.communicate()

        with open(filename, 'w') as f:
            f.write('# Mesh generated on: %s\n' % time)
            f.write('# Mesh generator version: %s\n' % commit)
            f.write(yaml.dump(self.original, default_flow_style=False))


    def sanity_check(self):
        sum_elems = self.n_te + self.n_back + self.n_front
        check_error(sum_elems % 4 == 0, "Number of radial elements must be a multiple of four")
        check_warn(sum_elems % 8 == 0, "Number of radial elements should be a multiple of eight")

        check_error(self.mesh_mode in {'2d', 'semi3d', '3d'},
                    "Invalid value for mesh_mode (valid: 2d, semi3d, 3d)")
        check_error(self.length_mode in {'extruded', 'uniform', 'double', 'triple'},
                    "Invalid value for length_mode (valid: extruded, uniform, double, triple)")
        check_error(self.order in {2, 3, 4}, "Order must be 2, 3 or 4")

        if self.mesh_mode != '2d':
            if self.length_mode == 'extruded':
                check_error(len(self.wingdef) == 1, "More than one wing definition in extruded mode")
                check_warn(self.n_base == 0, "n_base > 0 has no effect in extruded mode")
                check_warn(self.join_adds == 0, "join_adds > 0 has no effect in extruded mode")
                check_error(self.length > 0, "length must be positive in extruded mode")
            elif self.length_mode == 'uniform':
                check_error(len(self.wingdef) > 1, "Less than two wing definitions in uniform mode")
                check_warn(self.n_base == 0, "n_base > 0 has no effect in uniform mode")
            elif self.length_mode == 'double':
                check_error(len(self.wingdef) > 1, "Less than two wing definitions in double mode")
                check_error(self.n_length % 2 == 0, "n_length must be even in double mode")
                check_warn(self.n_base == 0, "n_base > 0 has no effect in double mode")
            elif self.length_mode == 'triple':
                check_error(len(self.wingdef) > 1, "Less than two wing definitions in triple mode")
                check_error(self.n_base > 0, "n_base should be positive in triple mode")
                check_error(self.n_length > self.n_base, "n_length <= n_base in triple mode")
                check_error((self.n_length - self.n_base) % 2 == 0,
                            "n_length - n_base should be even in triple mode")
        else:
            check_error(len(self.wingdef) == 1, "More than one wing definition in 2D mode")
            check_warn(self.join_adds == 0, "join_adds > 0 has no effect in extruded mode")

        for attr in ['sides', 'behind', 'ahead']:
            check_error(getattr(self, attr) >= 0, attr + ' < 0')
            check_error(getattr(self, attr) == 0 or
                        0 < getattr(self, 'p_' + attr) <= getattr(self, 'n_' + attr),
                        'Condition broken: 0 < p_' + attr + ' <= n_' + attr)
        check_error(0 < self.p_inner <= self.n_circle + self.n_square,
                    "Condition broken: 0 < p_inner < n_circle + n_square")

        for a, b in product(['in', 'out', 'hub', 'antihub', 'slip'], repeat=2):
            attr = a + '_' + b
            if hasattr(self, attr):
                check_error(getattr(self, attr) in [a, b],
                            "Invalid value for %s (valid: %s, %s)" % (attr, a, b))


    def radial_grading(self, length):
        fac = grading(length, self.len_bndlayer, self.n_bndlayers)
        return fac ** (1. / self.n_bndlayer)


    def angular_splits(self):
        n_quarter = (self.n_te + self.n_back + self.n_front) / 2
        n_small = n_quarter / 2
        n_large = n_quarter - n_small

        return [0, n_small, n_quarter, n_quarter + n_large, 2*n_quarter,
                2*n_quarter + n_small, 3*n_quarter, 3*n_quarter + n_large,
                4*n_quarter]


    def make_folder(self, folder):
        try:
            os.makedirs(folder)
        except OSError:
            pass


    def dump_g2files(self, folder, patches, always=False):
        if not (self.debug or always):
            return

        fn = 'out/{}'.format(folder)
        self.make_folder(fn)

        fn += '/{:03}.g2'
        for i, p in enumerate(patches):
            WriteG2(fn.format(i+1), p.objects())


    def dump_g2file(self, name, patches, always=False):
        if not (self.debug or always):
            return

        fn = 'out/{}.g2'.format(name)
        WriteG2(fn, patches)
