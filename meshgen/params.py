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
    'format': 'IFEM',
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
    print "Please study README.md for instructions."


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
        check_error(sum_elems % 4 == 0, "Number of angular elements must be a multiple of four")

        check_error(self.format in {'IFEM', 'OpenFOAM'},
                    "Invalid value for format (valid: IFEM, OpenFOAM)")
        check_error(self.mesh_mode in {'2d', 'semi3d', '3d'},
                    "Invalid value for mesh_mode (valid: 2d, semi3d, 3d)")
        check_error(self.length_mode in {'extruded', 'uniform', 'double', 'triple'},
                    "Invalid value for length_mode (valid: extruded, uniform, double, triple)")
        check_error(self.order in {2, 3, 4}, "Order must be 2, 3 or 4")

        if self.format == 'OpenFOAM':
            check_warn(self.mesh_mode != 'semi3d', "Semi3D output for OpenFOAM is incomplete (no beam)")
            check_error(self.order == 2, "OpenFOAM format requires order=2")
            check_error(not self.walldistance, "OpenFOAM format does not support wall distances")

        if self.mesh_mode != '2d':
            if self.length_mode == 'extruded':
                check_error(len(self.wingdef) == 1, "More than one wing definition in extruded mode")
                check_warn(self.n_base == 0, "n_base > 0 has no effect in extruded mode")
                check_warn(self.join_adds == 0, "join_adds > 0 has no effect in extruded mode")
                check_error(self.length > 0, "length must be positive in extruded mode")
            elif self.length_mode == 'uniform':
                check_error(len(self.wingdef) > 1, "Fewer than two wing definitions in uniform mode")
                check_warn(self.n_base == 0, "n_base > 0 has no effect in uniform mode")
            elif self.length_mode == 'double':
                check_error(len(self.wingdef) > 1, "Fewer than two wing definitions in double mode")
                check_error(self.n_length % 2 == 0, "n_length must be even in double mode")
                check_warn(self.n_base == 0, "n_base > 0 has no effect in double mode")
            elif self.length_mode == 'triple':
                check_error(len(self.wingdef) > 1, "Fewer than two wing definitions in triple mode")
                check_error(self.n_base > 0, "n_base should be positive in triple mode")
                check_error(self.n_length > self.n_base, "n_length <= n_base in triple mode")
                check_error((self.n_length - self.n_base) % 2 == 0,
                            "n_length - n_base should be even in triple mode")

            check_warn(self.n_length % self.p_length == 0,
                       "n_length should be a multiple of p_length for load balancing purposes")
        else:
            check_error(len(self.wingdef) == 1, "More than one wing definition in 2D mode")
            check_warn(self.join_adds == 0, "join_adds > 0 has no effect in extruded mode")

        for attr in ['sides', 'behind', 'ahead']:
            pn, nn = 'p_' + attr, 'n_' + attr
            check_error(getattr(self, attr) >= 0, attr + ' < 0')
            check_error(getattr(self, attr) == 0 or
                        0 < getattr(self, pn) <= getattr(self, nn),
                        "Condition broken: 0 < %s <= %s" % (pn, nn))
            check_warn(getattr(self, nn) % getattr(self, pn) == 0,
                       "%s should be a mulitple of %s for load balancing purposes" % (nn, pn))

        check_error(0 < self.p_inner <= self.n_circle + self.n_square,
                    "Condition broken: 0 < p_inner <= n_circle + n_square")
        check_warn((self.n_circle + self.n_square) % self.p_inner == 0,
                   "n_circle + n_square should be a multiple of p_inner for load balancing purposes")

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
