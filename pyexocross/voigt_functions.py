import numpy as np
from numba import jit

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

# def compute_fwhm_voigt_olivero(sigma, gamma):
#     from .constants import RT2LN2
#     import numpy as np
#     fwhm_g = 2*sigma*RT2LN2
#     fwhm_L = 2*gamma
#     return 0.5346*fwhm_L + np.sqrt(0.2166*fwhm_L**2 + fwhm_g**2)



# def whiting_method(v0, sigma, gamma):
#     fw_v = compute_fwhm_voigt_olivero(sigma, gamma)
#     fw_L = 2*gamma

#     v0


# def voigt_radis(v0, sigma, gamma):
#     wv = olivero_1977(2*sigma, 2*gamma)
#     wl = 2* gamma
#     return _whiting_jit(v0,wl, wv)*2.2

def fwhm(sigma, gamma):
    from .constants import RT2LN2
    fwhm_g = 2*sigma*RT2LN2
    fwhm_l = 2*gamma
    return fwhm_g, fwhm_l

def psuedo_voigt(v0, sigma, gamma):
    f_g, f_l = fwhm(sigma,gamma)

    f = (f_g**5 + 2.69269*f_g**4*f_l + 2.42843*f_g**3*f_l**2 +\
        4.47163*f_g**2*f_l**3 + 0.07842*f_g*f_l**4 + f_l**5)**0.2
    
    n = 1.36603*(f_l/f) -0.47717*(f_l/f)**2 + 0.11116*(f_l/f)**3

    return n*lorentz(v0,f/2) + (1-n)*gaussian(v0,f/2)

Voigt.add_voigt_function('psuedo', psuedo_voigt)

def gaussian(x, sigma):
    from .constants import SQRT2PI
    return np.exp(-x**2/(2*sigma**2))/(sigma*SQRT2PI)

def lorentz(x, gamma):
    return gamma/(np.pi*(x**2 + gamma**2))



def voigt_convolution(v0, sigma, gamma):
    from scipy.signal import fftconvolve
    a = gaussian(v0, sigma)
    b = lorentz(v0, gamma)
    return fftconvolve(a,b,mode='same')






    