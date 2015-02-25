import os
import xml.etree.ElementTree as xml

from collections import namedtuple
from math import sqrt
from shutil import rmtree

from GoTools import WriteG2
from GeoUtils.IO import ParseArgs

from utils import grading, gradspace


defaults = {
    # Input and output
    'wingfile': 'NREL_5_MW.xinp',
    'out': 'NREL',

    # Mesh generator
    'debug': False,
    'nprocs_mg': 6,

    # Global mesh parameters
    'order': 2,

    # Trailing edge
    'len_te': 2e-2,
    'len_te_cyl_fac': 50,

    # Angular resolution
    'n_te': 9,
    'n_back': 28,
    'n_front': 39,

    # Lengthwise resolution
    'length_mode': 'triple',
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
    'ahead': 0.2,
    'sides': 0.2,
    'n_behind': 36,
    'n_ahead': 2,
    'n_sides': 2,
}


def fix_floats(dct, keys=['z', 'theta', 'chord', 'ac', 'ao']):
    for k in keys:
        dct[k] = float(dct[k])
    return dct


Section = namedtuple('Section', ['z', 'theta', 'chord', 'ac', 'ao', 'foil'])


class Params(object):

    def __init__(self, args=[]):
        ParseArgs(args, defaults)

        wingdef = xml.parse(defaults['wingfile'])
        defaults['wingdef'] = [Section(**fix_floats(s.attrib))
                               for s in wingdef.getroot()]

        self.__dict__.update(defaults)

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


    def sanity_check(self):
        assert((self.n_te + self.n_back + self.n_front) % 4 == 0)

        assert(self.length_mode in ['extruded', 'uniform', 'double', 'triple'])

        if self.length_mode == 'extruded':
            assert(len(self.wingdef) == 1)
            assert(self.n_base == 0)
        elif self.length_mode == 'uniform':
            assert(len(self.wingdef) > 1)
            assert(self.n_base == 0)
        elif self.length_mode == 'double':
            assert(len(self.wingdef) > 1)
            assert(self.n_base == 0)
            assert(self.n_length % 2 == 0)
        elif self.length_mode == 'triple':
            assert(len(self.wingdef) > 1)
            assert(self.n_base > 0)
            assert(self.n_length > self.n_base)
            assert((self.n_length - self.n_base) % 2 == 0)
            assert(self.join_index > 0)


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
