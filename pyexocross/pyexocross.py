from .linelist import Linelist
from .voigt_functions import Voigt

class PyExocross:

    def __init__(self, linelist: Linelist, compute_voigt=True):
        self._linelist = linelist
        self._voigt = Voigt()
        self._compute_voigt = compute_voigt


    def set_voigt_function(self, voigt_function):
        self._voigt.set_voigt_function(voigt_function)
    
    def compute_xsec(self,wngrid,T,P,with_progress=True, chunksize=10000,
                     wing_cutoff=25.0, threshold=1e-34):
        import numpy as np
        from progress.counter import Counter
        import math
        bar = None
        _wngrid = np.sort(wngrid)

        total_transitions = self._linelist.totalTransitions
        num_elems = (int(math.ceil(total_transitions/chunksize))-1)*2
        if with_progress:
            bar=Counter('Processing transitions:')
        xsec = []

        out_grid = []
        if self._compute_voigt:
            xsec = np.zeros_like(wngrid)
        
        for v, I, gamma, doppler,count in self._linelist.transitions(_wngrid,
                                                               T, P,
                                                               wing_cutoff=wing_cutoff,
                                                               chunksize=chunksize,
                                                               threshold=threshold):
            if bar is not None:
                bar.next(n=count)
            if v is None or len(v) == 0:
                continue
            if self._compute_voigt:
                self._voigt.voigt(_wngrid, v, I, doppler, gamma,cutoff=wing_cutoff, out=xsec) 
            else:
                out_grid.append(v)
                xsec.append(I)
        
        if not self._compute_voigt:
            wngrid = np.concatenate(out_grid)
            xsec = np.concatenate(xsec)

        return wngrid, xsec