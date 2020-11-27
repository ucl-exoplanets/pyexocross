import numpy as np
from scipy.interpolate import interp1d
class ExomolPF:

    def __init__(self, filename):
        pfarray = np.loadtxt(filename)
        self._f = interp1d(pfarray[:,0], pfarray[:,1])
    def Q(self, temperature):
        return self._f(temperature)[()]



