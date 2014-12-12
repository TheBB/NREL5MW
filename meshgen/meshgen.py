from airfoils import AirFoil


class MeshGen(object):

    def __init__(self, params):
        self.params = params

    def load_airfoils(self):
        self.airfoils = [AirFoil(self.params, wd) for wd in self.params.wingdef]
        self.params.dump_g2files('airfoils_raw', self.airfoils)

    def resample_airfoils(self):
        for af in self.airfoils:
            af.resample(self.params)
        self.params.dump_g2files('airfoils_resampled', self.airfoils)
