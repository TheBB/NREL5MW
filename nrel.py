import sys

from meshgen.params import Params
from meshgen.meshgen import MeshGen

params = Params(sys.argv[1:])
gen = MeshGen(params)

gen.load_airfoils()
gen.resample_airfoils()
gen.resolve_join()
gen.resample_length()
gen.length_tfi()
