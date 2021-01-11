from pyexocross.exomol.exomollinelist import ExomolLinelist
from pyexocross.pyexocross import PyExocross
import numpy as np
from pyexocross.util import create_grid_res, convert_to_wavenumber
import matplotlib.pyplot as plt


if __name__=="__main__":

    T = np.linspace(100,500,5)
    P = np.logspace(-13,-3,5)

    wlgrid = create_grid_res(100000,0.5,1.0)[:,0]
    wngrid =convert_to_wavenumber(wlgrid[::-1],'um')

    #wngrid = np.linspace(100,10000,10000)

    hl = ExomolLinelist('/Users/ahmed/Documents/molecular_data/Exomol/TiO')
    #hl = ExomolLinelist('/Users/ahmed/Documents/Linelists/H2O')
    hl.add_default_broadener()


    pyexo = PyExocross(hl)

    t = 3000
    p = 3.0

    wn_approx,xsec_approx = pyexo.compute_xsec_parallel(wngrid,t,p, chunksize=1000000, threshold=0.0, wing_cutoff=25.0, max_workers=2,max_jobs=4)
    print('\n SCIPY done\n')
    pyexo.set_voigt_function('psuedo')
    wn_acc,xsec_acc = pyexo.compute_xsec_parallel(wngrid,t,p, chunksize=1000000, threshold=0.0, wing_cutoff=25.0, max_workers=2,max_jobs=40)

    #wn,xsec = pyexo.compute_xsec(wngrid,t,p, chunksize=1000, threshold=0.0, wing_cutoff=25.0)
    plt.figure()
    plt.plot(wn_approx,xsec_approx,label='scipy')
    plt.plot(wn_acc,xsec_acc,label='psuedo')
    plt.xlabel(r'Wavenumber cm$^{-1}$')
    plt.ylabel(r'Cross-section cm$^{2}$/molecule')
    plt.yscale('log')
    plt.legend()
    plt.show()

