import astropy.units as u

def convert_to_wavenumber(value, units):
    from_unit = u.Unit(units)
    conversion = (value * from_unit).to(u.k, equivalencies=u.spectral())

    return conversion.value


