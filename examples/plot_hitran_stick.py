from pyexocross.hitran.hitran import HITRANLinelist
from pyexocross.pyexocross import PyExocross
import numpy as np
from pyexocross.util import create_grid_res, convert_to_wavenumber
from pyexocross.writer.hdf5writer import HDF5Writer
import matplotlib.pyplot as plt
T = np.linspace(100,500,5)
P = np.logspace(-13,-3,5)


wngrid = np.array([0.0,10000])

hl = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CO2/5fc4ccb6.par')
#hl = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CO2/CO2all.par')
hl.add_self_broadener(ratio=0.99)
hl.add_air_broadener(ratio=0.01)


pyexo = PyExocross(hl,compute_voigt=False)

t = 296
p = 1
plt.figure()
wn,xsec = pyexo.compute_xsec(wngrid,t,p, chunksize=1000, threshold=1e-25)
print('\n',len(xsec))
plt.stem(10000/wn,xsec,linefmt='-',markerfmt=" ")
plt.xlabel(r'Wavenumber cm$^{-1}$')
plt.ylabel(r'Absolute Intensity cm/molecule')
plt.yscale('log')
plt.xscale('log')
plt.show()

