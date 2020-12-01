from pyexocross.hitran.hitran import HITRANLinelist
from pyexocross.pyexocross import PyExocross
import numpy as np
from pyexocross.util import create_grid_res, convert_to_wavenumber
from pyexocross.writer.hdf5writer import HDF5Writer

nT = 5
nP = 15

T = np.linspace(150,300,nT)
P = np.logspace(-14,1,nP)
wlgrid = create_grid_res(100000,1.0,10)[:,0]
wngrid =convert_to_wavenumber(wlgrid[::-1],'um')

arr = np.loadtxt('/Users/ahmed/Documents/repos/pyexocross/examples/psg1.dat')
wngrid = np.sort(arr[:,0])

hl = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CO2/12C16O2.par')
hl.add_self_broadener(ratio=0.96)
hl.add_air_broadener(ratio=0.04)
print(wngrid)
pyexo = PyExocross(hl, compute_voigt=True)

total = nP*nT
count = 0
with HDF5Writer('xsecs/all_ab.h5','CO2',T,P) as f:
    for t in T:
        for p in P:
            print(f'\n\nT={t}  P={p} [{count}/{total}]')
            wn,xsec = pyexo.compute_xsec(wngrid,t,p, chunksize=30000, threshold=0.0, wing_cutoff=25.0)
            f.add_cross_section(t,p,wn,xsec)
            count += 1


