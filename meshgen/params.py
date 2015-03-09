# -*- coding: utf-8 -*-

import os
import sys
import xml.etree.ElementTree as xml

from itertools import product, combinations, permutations, chain
from collections import namedtuple
from datetime import datetime
from math import sqrt
from shutil import rmtree
import subprocess
import yaml

from GoTools import WriteG2, SetProcessorCount
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

    # Flow and fluid characteristics
    'Re': 10000.0,
    'rho': 1.0,
    'velocity': 1.0,

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
    'p_rigid': 0,
}


def usage():
    print "Please study README.md for instructions."


def fix_floats(dct, keys=['z', 'theta', 'chord', 'ac', 'ao']):
    """Coerces the given keys in the dictionary to floats."""
    for k in keys:
        dct[k] = float(dct[k])
    return dct


class Colors(object):
    """Enum of colors used for ANSI terminal output."""
    WARNING = '\033[93m\033[1m'
    ERROR = '\033[91m\033[1m'
    END = '\033[0m'


def check_warn(test, msg):
    """Produces a warning message if the test fails."""
    if not test:
        print Colors.WARNING + "WARNING: " + Colors.END + msg


def check_error(test, msg):
    """Produces an error message and exits if the test fails."""
    if not test:
        print Colors.ERROR + "ERROR: " + Colors.END + msg
        sys.exit(0)


Section = namedtuple('Section', ['z', 'theta', 'chord', 'ac', 'ao', 'foil'])


class Params(object):
    """This class holds all the parameters."""

    def __init__(self, args=[]):
        # Parse all the arguments using the defaults
        ParseArgs(args, defaults)

        # Store the parsed arguments in their original state (no postprocessing)
        self.original = defaults

        # Update the namespace
        self.__dict__.update(defaults)

        # Parse the wing definition file
        wingdef = xml.parse(self.wingfile)
        self.wingdef = [Section(**fix_floats(s.attrib)) for s in wingdef.getroot()]

        # Perform some minor postprocessing
        self._postprocess()
        self._easy_extensions()

        # Warnings and errors in case something is wrong
        self._sanity_check()

        # Act on the parameters if necessary
        self._act()

        # Write a status message
        s = {'blade': 'Blade mode',
             '2d': '2D mode',
             'semi3d': 'Semi3D mode (%i planes)' % (self.n_length + 1),
             '3d': '3D mode'}[self.mesh_mode]
        s += ' -- ' + {2: 'linear', 3: 'quadratic', 4: 'cubic'}[self.order] + ' geometry'
        s += ' -- ' + self.format + ' format'
        print s

    def out_yaml(self):
        """Writes a YAML file with the original parameters."""
        filename = os.path.join(self.out, 'parameters.yaml')

        time = datetime.now().isoformat()
        proc = subprocess.Popen(['git', '-C', os.path.abspath(os.path.dirname(__file__)),
                                 'rev-parse', 'HEAD'], stdout=subprocess.PIPE)
        commit, _ = proc.communicate()

        with open(filename, 'w') as f:
            f.write('# Mesh generated on: %s\n' % time)
            f.write('# Mesh generator version: %s\n' % commit)
            f.write(yaml.dump(self.original, default_flow_style=False))

    def radial_grading(self, length):
        """Computes the radial grading factor needed for the given length to the circle."""
        fac = grading(length, self.len_bndlayer, self.n_bndlayers)
        return fac ** (1. / self.n_bndlayer)

    def angular_splits(self):
        """Computes the sizes of each angular part."""
        n_quarter = (self.n_te + self.n_back + self.n_front) / 2
        n_small = n_quarter / 2
        n_large = n_quarter - n_small

        return [0, n_small, n_quarter, n_quarter + n_large, 2*n_quarter,
                2*n_quarter + n_small, 3*n_quarter, 3*n_quarter + n_large,
                4*n_quarter]

    def dump_g2files(self, folder, patches, always=False):
        """Dump debug g2-files. The objects in the list patches must support the .objects()
        method. This method only does something if in debug mode, or if always is True."""
        if not (self.debug or always):
            return

        fn = 'out/{:02}_{}'.format(self._num_out, folder)
        self._make_folder(fn)

        fn += '/{:03}.g2'
        for i, p in enumerate(patches):
            WriteG2(fn.format(i+1), list(p.objects()))

        self._num_out += 1

    def dump_g2file(self, name, patches, always=False):
        """Dump a single debug file with the given patches. This method only does something if in
        debug mode, or if always is True."""
        if not (self.debug or always):
            return

        fn = 'out/{:02}_{}.g2'.format(self._num_out, name)
        WriteG2(fn, patches)

        self._num_out += 1

    def postprocess_xinp(self, path):
        """Performs postprocessing on an IFEM xinp file."""
        # Output flow characteristics
        mu = self.rho * self.velocity * self.len_char / self.Re

        with open(path, 'a') as f:
            f.write('\n<stokes>\n')
            f.write('  <fluidproperties mu="%e" rho="%e" />\n' % (mu, self.rho))
            f.write('</stokes>\n')

    def _postprocess(self):
        """Performs all the postprocessing."""
        self.len_te_cyl = self.len_te_cyl_fac * self.len_te
        self.len_char = min(wd.chord for wd in self.wingdef)
        self.len_bndlayer = self.len_char / sqrt(self.Re)
        self.n_bndlayers = float(self.n_circle) / self.n_bndlayer

        self.behind *= self.radius
        self.ahead *= self.radius
        self.sides *= self.radius

        self._num_out = 1

    def _easy_extensions(self):
        """This method produces ext_xyz_not_abc-names for all subsets xyz and abc of possible
        extensions.  This is basically just for simplicity. E.g.
        p.behind > 0 => p.ext_b
        p.behind > 0 and p.sides > 0 => p.ext_bs, p.ext_sb
        p.behind > 0 and p.sides == 0 => p.ext_b_not_s
        p.behind == 0 and p.sides == 0 => p.ext_not_bs, p.ext_not_sb
        etc. Because I can."""
        def subsets(s):
            for length in xrange(0, len(s) + 1):
                for comb in combinations(s, length):
                    yield comb

        extensions = {'behind', 'ahead', 'sides'}
        for on_set in subsets(extensions):
            remaining = extensions - set(on_set)
            for off_set in subsets(remaining):
                if not on_set and not off_set:
                    continue
                for on, off in product(permutations(on_set), permutations(off_set)):
                    name = 'ext'
                    if on:
                        name += '_' + ''.join(attr[0] for attr in on)
                    if off:
                        name += '_not_' + ''.join(attr[0] for attr in off)
                    val = all(getattr(self, attr) > 0 for attr in on)
                    val = val and all(getattr(self, attr) == 0 for attr in off)
                    setattr(self, name, val)

    def _act(self):
        """Performs all the pre-actions necessary depending on the parameters given."""
        SetProcessorCount(self.nprocs_mg)

        if self.debug:
            rmtree('out', ignore_errors=True)
            self._make_folder('out')

        rmtree(self.out, ignore_errors=True)
        self._make_folder(self.out)

    def _sanity_check(self):
        """Performs a basic sanity check of the parameters."""
        # Check parameters that must be in a set of allowed values
        check_error(self.format in {'IFEM', 'OpenFOAM'},
                    "Invalid value for format (valid: IFEM, OpenFOAM)")
        check_error(self.mesh_mode in {'blade', '2d', 'semi3d', '3d'},
                    "Invalid value for mesh_mode (valid: blade, 2d, semi3d, 3d)")
        check_error(self.length_mode in {'extruded', 'uniform', 'double', 'triple'},
                    "Invalid value for length_mode (valid: extruded, uniform, double, triple)")
        check_error(self.order in {2, 3, 4}, "Order must be 2, 3 or 4")

        for a, b in product(['in', 'out', 'hub', 'antihub', 'slip'], repeat=2):
            attr = a + '_' + b
            if hasattr(self, attr):
                check_error(getattr(self, attr) in [a, b],
                            "Invalid value for %s (valid: %s, %s)" % (attr, a, b))
        # Output format checks
        if self.format == 'OpenFOAM':
            check_warn(self.mesh_mode != 'semi3d', "Semi3D output for OpenFOAM is incomplete (no beam)")
            check_error(self.order == 2, "OpenFOAM format requires order=2")
            check_error(not self.walldistance, "OpenFOAM format does not support wall distances")

        # Checks for mesh_mode and length_mode
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

        # Checks for extensions
        for attr in ['sides', 'behind', 'ahead']:
            pn, nn = 'p_' + attr, 'n_' + attr
            check_error(getattr(self, attr) >= 0, attr + ' < 0')
            check_error(getattr(self, attr) == 0 or
                        0 < getattr(self, pn) <= getattr(self, nn),
                        "Condition broken: 0 < %s <= %s" % (pn, nn))
            check_warn(getattr(self, nn) % getattr(self, pn) == 0,
                       "%s should be a mulitple of %s for load balancing purposes" % (nn, pn))

        # Angular refinement checks
        sum_elems = self.n_te + self.n_back + self.n_front
        check_error(sum_elems % 4 == 0, "Number of angular elements must be a multiple of four")

        # Radial refinement checks
        check_error(0 < self.p_inner <= self.n_circle + self.n_square,
                    "Condition broken: 0 < p_inner <= n_circle + n_square")
        check_warn((self.n_circle + self.n_square) % self.p_inner == 0,
                   "n_circle + n_square should be a multiple of p_inner for load balancing purposes")

        # Rigid part
        check_error(0 <= self.p_rigid <= self.p_inner, "Condition broken: 0 <= p_rigid <= p_inner")

    def _make_folder(self, folder):
        """Make a folder if it doesn't exist."""
        try:
            os.makedirs(folder)
        except OSError:
            pass
