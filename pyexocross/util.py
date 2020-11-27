import astropy.units as u

def convert_to_wavenumber(value, units):
    from_unit = u.Unit(units)
    conversion = (value * from_unit).to(u.k, equivalencies=u.spectral())

    return conversion.value


def compute_gamma(gamma0,n0, T, P, t0, p0):
    
    return gamma0*((t0/T)**n0)*(P/p0)