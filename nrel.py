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

gen.subdivide_airfoils()
gen.lower_order()

if params.length_mode == 'extruded':
    gen.extrude()

if params.mesh_mode == 'semi3d':
    gen.output_semi3d()
elif params.mesh_mode == '2d':
    gen.output_2d()
