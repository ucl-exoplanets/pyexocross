from pyexocross.hitran.hitran import HITRANLinelist
from pyexocross.pyexocross import PyExocross
from pyexocross.exomol.exomolbroads import ExomolBroadener
import numpy as np
from pyexocross.util import create_grid_res, convert_to_wavenumber
from pyexocross.writer.hdf5writer import HDF5Writer
import matplotlib.pyplot as plt

wngrid = 10000/create_grid_res(15000,1.1,2.0)[::-1,0]

#hl_h2o = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/H2O/H2O.par')
hl_h2o= HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CH4/CH4.par')
#hl = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CO2/12C16O2.par')
h2_h2o = ExomolBroadener(0.0209,0.027,filename='/Users/ahmed/Documents/molecular_data/HITRAN/CH4/1H2-16O__H2.broad',species='H2')
he_h2o = ExomolBroadener(0.0042,0.20,filename='/Users/ahmed/Documents/molecular_data/HITRAN/CH4/1H2-16O__He.broad',species='He')

hl_h2o.add_broadener(h2_h2o,ratio=0.704)
hl_h2o.add_broadener(he_h2o,ratio=0.121)
hl_h2o.add_self_broadener(ratio=0.1)

# h2_ch4 = ExomolBroadener(0.0603,0.5,filename='/Users/ahmed/Documents/molecular_data/HITRAN/CH4/12C-1H4__H2.broad',species='H2')
# he_ch4 = ExomolBroadener(0.0382,0.30,filename='/Users/ahmed/Documents/molecular_data/HITRAN/CH4/12C-1H4__He.broad',species='He')

# hl_ch4.add_broadener(h2_ch4,ratio=0.83)
# hl_ch4.add_broadener(he_ch4,ratio=0.17)

pyexo_h2o = PyExocross(hl_h2o)
#pyexo_ch4 = PyExocross(hl_ch4)


t = 200
p = 1.0
if __name__ == "__main__":
    wn_h2o_self,xsec_h2o_self = pyexo_h2o.compute_xsec_parallel(wngrid,t,p, chunksize=1000, threshold=0.0, wing_cutoff=25.0,max_workers=2)
    hl_h2o.set_broadener_ratio('self',ratio=1e-10)
    wn_h2o,xsec_h2o = pyexo_h2o.compute_xsec_parallel(wngrid,t,p, chunksize=1000, threshold=0.0, wing_cutoff=25.0,max_workers=2)
    #wn_ch4,xsec_ch4 = pyexo_ch4.compute_xsec(wngrid,t,p, chunksize=100, threshold=0.0, wing_cutoff=25.0)
    plt.figure()
    # plt.plot(wn,xsec,label='pyexo')
    plt.plot(wn_h2o_self,xsec_h2o_self,label='H2O self')
    plt.plot(wn_h2o,xsec_h2o,label='H2O')
    #plt.plot(10000/wn_ch4,xsec_ch4,label='CH4')
    plt.xlabel(r'Wavelength um')
    plt.ylabel(r'Cross-section cm$^{2}$/molecule')
    plt.yscale('log')
    plt.legend()
    plt.show()

