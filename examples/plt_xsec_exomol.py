from pyexocross.exomol.exomollinelist import ExomolLinelist
from pyexocross.pyexocross import PyExocross
import numpy as np
from pyexocross.util import create_grid_res, convert_to_wavenumber
import matplotlib.pyplot as plt
T = np.linspace(100,500,5)
P = np.logspace(-13,-3,5)

#wlgrid = create_grid_res(1000,0.1,10)[:,0]
#wngrid =convert_to_wavenumber(wlgrid[::-1],'um')

wngrid = np.linspace(2000,3000,10000)

hl = ExomolLinelist('/Users/ahmed/Documents/molecular_data/Exomol/TiO')
hl.add_default_broadener()


pyexo = PyExocross(hl)

t = 296
p = 5
plt.figure()
wn,xsec = pyexo.compute_xsec(wngrid,t,p, chunksize=100000,threshold=1e-28)
plt.plot(wn,xsec)
plt.xlabel(r'Wavenumber cm$^{-1}$')
plt.ylabel(r'Cross-section cm$^{2}$/molecule')
plt.show()

