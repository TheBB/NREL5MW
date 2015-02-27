import sys

from meshgen.params import Params, usage
from meshgen.meshgen import MeshGen

if 'help' in sys.argv:
    usage()
    sys.exit(0)

params = Params(sys.argv[1:])
gen = MeshGen(params)

gen.load_airfoils()
gen.resample_airfoils()
gen.resolve_join()
gen.resample_length()
gen.fill_airfoils()
gen.fill_sides()

if params.mesh_mode == 'semi3d':
    gen.subdivide_airfoils()
    gen.lower_order()
    gen.output_planes()

