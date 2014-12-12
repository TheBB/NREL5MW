from airfoils import AirFoil


class MeshGen(object):

    def __init__(self, params):
        self.params = params


    def load_airfoils(self):
        self.airfoils = [AirFoil(params=self.params, wd=wd) for wd in self.params.wingdef]
        self.params.dump_g2files('airfoils_raw', self.airfoils)


    def resample_airfoils(self):
        for af in self.airfoils:
            af.resample(self.params.n_te, self.params.n_back, self.params.n_front)
        self.params.dump_g2files('airfoils_resampled', self.airfoils)


    def resolve_join(self):
        n = self.params.join_adds
        idx = self.params.join_index

        if n == 0:
            return

        for i in xrange(n):
            new = AirFoil.average(self.airfoils[idx], self.airfoils[idx+1])
            self.airfoils.insert(idx + 1, new)

        for i in xrange(n):
            new = AirFoil.average(self.airfoils[idx+i-1], self.airfoils[idx+i])
            self.airfoils.insert(idx + i, new)

        del self.airfoils[idx + n]

        self.params.dump_g2files('airfoils_joined', self.airfoils)
