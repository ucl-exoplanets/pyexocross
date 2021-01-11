from pyexocross.hitran.hitran import HITRANLinelist
from pyexocross.pyexocross import PyExocross
from pyexocross.exomol.exomolbroads import ExomolBroadener
import numpy as np
from pyexocross.util import create_grid_res, convert_to_wavenumber
from pyexocross.writer.hdf5writer import HDF5Writer

if __name__ =="__main__":
    nP = 22

    T = np.array([200,300])
    nT = len(T)
    P = np.logspace(-5,2,22)
    wlgrid = create_grid_res(1000000,1.3,1.5)[:,0]
    wngrid =convert_to_wavenumber(wlgrid[::-1],'um')
    wngrid = np.sort(wngrid)

    hl_h2o = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/H2O/H2O.par')
    #hl_ch4 = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CH4/CH4.par')
    #hl = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CO2/12C16O2.par')
    h2_h2o = ExomolBroadener(0.0209,0.027,filename='/Users/ahmed/Documents/molecular_data/HITRAN/H2O/1H2-16O__H2.broad',species='H2')
    he_h2o = ExomolBroadener(0.0042,0.20,filename='/Users/ahmed/Documents/molecular_data/HITRAN/H2O/1H2-16O__He.broad',species='He')

    hl_h2o.add_broadener(h2_h2o,ratio=0.704)
    hl_h2o.add_broadener(he_h2o,ratio=0.121)
    hl_h2o.add_self_broadener(ratio=0.107)
    pyexo_h2o = PyExocross(hl_h2o, compute_voigt=True)

    total = nP*nT
    count = 0
    with HDF5Writer('xsecs/H2O_self_hires.h5','H2O',T,P) as f:
        for t in T:
            for p in P:
                print(f'\n\nT={t}  P={p} [{count}/{total}]')
                wn,xsec = pyexo_h2o.compute_xsec_parallel(wngrid,t,p, chunksize=1000, threshold=0.0, wing_cutoff=25.0,max_workers=2)
                f.add_cross_section(t,p,wn,xsec)
                count += 1


