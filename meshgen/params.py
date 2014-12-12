import os
import xml.etree.ElementTree as xml

from collections import namedtuple
from shutil import rmtree

from GoTools import WriteG2
from GeoUtils.IO import ParseArgs


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
}


def fix_floats(dct, keys=['z', 'theta', 'chord', 'ac', 'ao']):
    for k in keys:
        dct[k] = float(dct[k])
    return dct


class Params(object):

    def __init__(self, args=[]):
        ParseArgs(args, defaults)

        wingdefs = xml.parse(defaults['wingfile'])
        defaults['wingdef'] = [
            namedtuple('Section', s.attrib.keys())(**fix_floats(s.attrib))
            for s in wingdefs.getroot()
        ]

        for key, value in defaults.iteritems():
            setattr(self, key, value)

        self.len_te_cyl = self.len_te_cyl_fac * self.len_te

        self.sanity_check()

        rmtree('out', ignore_errors=True)
        self.make_folder('out')


    def sanity_check(self):
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
