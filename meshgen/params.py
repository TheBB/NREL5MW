import os
import sys
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
    'nprocs_mg': 4,

    # Global mesh parameters
    'order': 2,

    # Trailing edge
    'len_te': 2e-2,
    'len_te_cyl_fac': 50,
}


def fix_floats(dct, keys=['z', 'theta', 'chord', 'ac', 'ao']):
    for k in keys:
        dct[k] = float(dct[k])
    return dct


class Params(object):

    def __init__(self):
        ParseArgs(sys.argv[1:], defaults)

        wingdefs = xml.parse(defaults['wingfile'])
        defaults['wingdef'] = [
            namedtuple('Section', s.attrib.keys())(**fix_floats(s.attrib))
            for s in wingdefs.getroot()
        ]

        for key, value in defaults.iteritems():
            setattr(self, key, value)

        self.len_te_cyl = self.len_te_cyl_fac * self.len_te

        rmtree('out', ignore_errors=True)
        self.make_folder('out')

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
