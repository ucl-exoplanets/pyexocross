from pyexocross.hitran.hitran import HITRANLinelist
from pyexocross.pyexocross import PyExocross
from pyexocross.exomol.exomolbroads import ExomolBroadener
import numpy as np
from pyexocross.util import create_grid_res, convert_to_wavenumber
from pyexocross.writer.hdf5writer import HDF5Writer
import matplotlib.pyplot as plt


if __name__ == "__main__":
    arr = np.loadtxt('/Users/ahmed/Documents/repos/pyexocross/examples/PSGtest.dat')
    #arr = np.loadtxt('/Users/ahmed/Documents/repos/pyexocross/examples/psg1.dat')
    #arr = np.loadtxt('/Users/ahmed/Documents/repos/pyexocross/examples/psglowpress.dat')
    #arr = np.loadtxt('/Users/ahmed/Documents/repos/pyexocross/examples/hightemplowpress.dat')
    #wngrid = np.sort(arr[:,0])
    wngrid=np.linspace(0,10000, 1000)
    spec = arr[:,1]
    hl = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CO2/CO2all.par')
    #hl = HITRANLinelist('/Users/ahmed/Documents/molecular_data/HITRAN/CO2/12C16O2.par')
    hl.add_self_broadener(ratio=1.0)


    pyexo = PyExocross(hl)


    t = 296
    p = 1

    wn_approx,xsec_approx = pyexo.compute_xsec_parallel(wngrid,t,p, chunksize=1000, threshold=0.0, wing_cutoff=25.0, max_workers=2,max_jobs=40)
    print('\n SCIPY done\n')
    pyexo.set_voigt_function('psuedo')
    wn_acc,xsec_acc = pyexo.compute_xsec_parallel(wngrid,t,p, chunksize=1000, threshold=0.0, wing_cutoff=25.0, max_workers=2,max_jobs=40)

    #wn,xsec = pyexo.compute_xsec(wngrid,t,p, chunksize=1000, threshold=0.0, wing_cutoff=25.0)
    plt.figure()
    plt.plot(wn_approx,xsec_approx,label='scipy')
    plt.plot(wn_acc,xsec_acc,label='psuedo')
    #plt.plot(arr[:,0], arr[:,1],label='psg')
    plt.xlabel(r'Wavenumber cm$^{-1}$')
    plt.ylabel(r'Cross-section cm$^{2}$/molecule')
    plt.yscale('log')
    plt.legend()
    plt.show()

