from pyexocross.hitran.hitran import HITRANLinelist
from pyexocross.pyexocross import PyExocross
import numpy as np
from pyexocross.util import create_grid_res, convert_to_wavenumber
from pyexocross.writer.hdf5writer import HDF5Writer
import matplotlib.pyplot as plt

arr = np.loadtxt('/Users/ahmed/Documents/repos/pyexocross/examples/PSGtest.dat')
#arr = np.loadtxt('/Users/ahmed/Documents/repos/pyexocross/examples/psg1.dat')
#arr = np.loadtxt('/Users/ahmed/Documents/repos/pyexocross/examples/psglowpress.dat')
#arr = np.loadtxt('/Users/ahmed/Documents/repos/pyexocross/examples/hightemplowpress.dat')
wngrid = np.sort(arr[:,0])
#wngrid=np.arange(3200,3800, 0.30382)
spec = arr[:,1]
hl = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CO2/CO2all.par')
#hl = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CO2/12C16O2.par')
hl.add_self_broadener(ratio=1.0)


pyexo = PyExocross(hl)


t = 296
p = 1

wn,xsec = pyexo.compute_xsec(wngrid,t,p, chunksize=1000, threshold=0.0, wing_cutoff=25.0)
plt.figure()
plt.plot(wn,xsec,label='pyexo')
plt.plot(arr[:,0], arr[:,1],label='psg')
plt.xlabel(r'Wavenumber cm$^{-1}$')
plt.ylabel(r'Cross-section cm$^{2}$/molecule')
plt.yscale('log')
plt.legend()
plt.show()

