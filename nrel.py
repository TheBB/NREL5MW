import sys

from meshgen.params import Params
from meshgen.meshgen import MeshGen

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
    gen.output_planes()
