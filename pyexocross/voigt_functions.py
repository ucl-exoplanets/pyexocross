import numpy as np


class Voigt:

    available_voigt ={
    }

    @classmethod
    def add_voigt_function(cls,name, func):
        if name not in cls.available_voigt:
            cls.available_voigt[name] =func

    def __init__(self, voigt_function='scipy', accurate_sum=False):
        self.set_voigt_function(voigt_function)
        self.use_accurate_sum(accurate_sum)

    def set_voigt_function(self, voigt_function):
        self._f = self.available_voigt[voigt_function]

    def use_accurate_sum(self, accurate_sum):
        self._sum = np.sum
        if accurate_sum:
            from accupy import fsum
            self._sum = fsum
    def voigt(self, wngrid, v, I, doppler, lorentz, cutoff=25.0, out=None):
        from .constants import RT2LN2
        min_v = v.min()-cutoff
        max_v = v.max()+cutoff

        start = max(0,wngrid.searchsorted(v.min()-cutoff)-1)
        end = min(wngrid.searchsorted(v.max()+cutoff),len(wngrid))
        res = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            x = wngrid[i] - v
            fil = np.abs(x) <=cutoff
            x=x[fil]
            sigma = doppler[fil]/RT2LN2
            gamma = lorentz[fil]
            res[idx] = self._sum(self._f(x,sigma,gamma)*I[fil])



        # x = wngrid[start:end,None] - v[None,:]

        # sigma = doppler/RT2LN2
        # gamma = lorentz
        # res = np.sum(self._f(x,sigma[None,:],gamma[None,:])*I[None,:],axis=1)
        if out is None:
            return res, start, end
        else:
            out[start:end] += res

def voigt_scipy(v0, sigma, gamma):
    from scipy.special import voigt_profile
    return voigt_profile(v0, sigma, gamma)
Voigt.add_voigt_function('scipy', voigt_scipy)

        









    