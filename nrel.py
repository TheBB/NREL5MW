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

if params.length_mode != 'extruded' and params.mesh_mode != '2d':
    gen.resample_length()

if params.mesh_mode == 'blade':
    if params.length_mode == 'extruded':
        gen.extrude()
    gen.loft_blade()
else:
    gen.make_slices()
    gen.extend()

    gen.subdivide_slices()

    if params.length_mode == 'extruded':
        gen.extrude()

    if params.mesh_mode == '3d':
        gen.loft_slices()

gen.lower_order()
gen.output()
