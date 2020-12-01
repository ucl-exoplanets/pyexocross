import astropy.units as u
import numpy as np
def convert_to_wavenumber(value, units):
    from_unit = u.Unit(units)
    conversion = (value * from_unit).to(u.k, equivalencies=u.spectral())

    return conversion.value


def compute_gamma(gamma0,n0, T, P, t0, p0):
    
    return gamma0*((t0/T)**n0)*(P/p0)

def create_grid_res(resolution, wave_min, wave_max):
    #
    # R = l/Dl
    # l = (l-1)+Dl/2 + (Dl-1)/2
    #
    # --> (R - 1/2)*Dl = (l-1) + (Dl-1)/2
    #
    # 
    wave_list = []
    width_list = []
    wave = wave_min
    width = wave/resolution    
    
    while wave < wave_max:
        width = wave / (resolution - 0.5) + width/2/(resolution - 0.5)
        wave = resolution * width 
        width_list.append(width)
        wave_list.append(wave)

    return np.array((wave_list ,width_list)).T