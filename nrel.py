from meshgen.params import Params
from meshgen.meshgen import MeshGen

params = Params()
gen = MeshGen(params)

gen.load_airfoils()
gen.resample_airfoils()
