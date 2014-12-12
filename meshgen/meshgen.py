from airfoils import AirFoil


class MeshGen(object):

    def __init__(self, params):
        self.params = params

    def load_airfoils(self):
        self.airfoils = [AirFoil(self.params, wd) for wd in self.params.wingdef]
        self.params.dump_g2files('airfoils', self.airfoils)
