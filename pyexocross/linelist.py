
class Linelist:

    @property
    def totalTransitions(self):
        raise NotImplementedError

    def transitions(self, wngrid, temperature, pressure, wing_cutoff=25.0, chunksize=10000, threshold=1e-34):
        raise NotImplementedError


