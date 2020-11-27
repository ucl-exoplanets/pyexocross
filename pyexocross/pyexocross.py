from .linelist import Linelist
from .voigt_functions import Voigt

class PyExocross:

    def __init__(self, linelist: Linelist):
        self._linelist = linelist
        self._voigt = Voigt()


    def set_voigt_function(self, voigt_function):
        self._voigt.set_voigt_function(voigt_function)
    
    def compute_xsec(self,wngrid,T,P,with_progress=True, chunksize=10000,
                     wing_cutoff=25.0, threshold=1e-34):
        import numpy as np
        from progress.bar import Bar
        import math
        bar = None
        total_transitions = self._linelist.totalTransitions
        num_elems = int(math.ceil(total_transitions/chunksize))
        if with_progress:
            bar = Bar('Computing', max=num_elems)
        
        xsec = np.zeros_like(wngrid)
        wngrid = wngrid
        for v, I, gamma, doppler in self._linelist.transitions(wngrid,
                                                               T, P,
                                                               wing_cutoff=25.0,
                                                               chunksize=chunksize,
                                                               threshold=0.0):
            if bar is not None:
                bar.next()
            if v is None or len(v) == 0:
                continue
            self._voigt.voigt(wngrid, v, I, doppler, gamma, out=xsec) 


        return wngrid, xsec