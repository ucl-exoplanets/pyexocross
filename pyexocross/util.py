import astropy.units as u
import numpy as np
def convert_to_wavenumber(value, units):
    from_unit = u.Unit(units)
    conversion = (value * from_unit).to(u.k, equivalencies=u.spectral())

    return conversion.value


def compute_gamma(gamma0,n0, T, P, t0, p0):
    
    return gamma0*((t0/T)**n0)*(P/p0)


def sanitize_molecule_string(molecule):
    import re
    """
    Cleans a molecule string to match up
    with molecule naming in TauREx3.

    e.g:

    H2O -> H2O

    1H2-16O -> H2O

    Parameters
    ----------
    molecule: str
        Molecule to sanitize
    
    Returns
    -------
    str:
        Sanitized name

    """
    return ''.join([''.join(s) for s in
                    re.findall('([A-Z][a-z]?)([0-9]*)', molecule)])

def conversion_factor(from_unit, to_unit):
    import astropy.units as u

    try:
        from_conv = u.Unit(from_unit)
    except:
        from_conv = u.Unit(from_unit, format="cds")

    try:
        to_conv = u.Unit(to_unit)
    except:
        to_conv = u.Unit(to_unit, format="cds")

    return from_conv.to(to_unit)

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


def doppler_broad(freq, mass, temperature):
    from .constants import KBOLTZ, SPDLIGT, AVGNO
    import math

    return np.sqrt(2*KBOLTZ*math.log(2)*AVGNO*temperature/mass)*freq/SPDLIGT