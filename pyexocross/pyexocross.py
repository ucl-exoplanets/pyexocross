from .linelist import Linelist
from .voigt_functions import Voigt

def parallel_voigt(args, wing_cutoff=25.0, wngrid=None):
    v, I, gamma, doppler,count = args
    if v is None:
        return None,None,None,None
    voigt = Voigt()
    return *voigt.voigt(wngrid, v, I, doppler, gamma,cutoff=wing_cutoff),count 


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
        
        import math
        bar = None
        _wngrid = np.sort(wngrid)

        total_transitions = self._linelist.totalTransitions
        num_elems = (int(math.ceil(total_transitions/chunksize))-1)*2
        itera = self._linelist.transitions(_wngrid,T, P,wing_cutoff=wing_cutoff,
                                                               chunksize=chunksize,
                                                               threshold=threshold)

        xsec = []

        out_grid = []
        if self._compute_voigt:
            xsec = np.zeros_like(wngrid)
        from tqdm import tqdm

        with tqdm() as t:
            for v, I, gamma, doppler,count in itera:

                if v is None or len(v) == 0:
                    continue
                t.update(count)
                if self._compute_voigt:
                    self._voigt.voigt(_wngrid, v, I, doppler, gamma,cutoff=wing_cutoff, out=xsec) 
                else:
                    out_grid.append(v)
                    xsec.append(I)
            
            if not self._compute_voigt:
                wngrid = np.concatenate(out_grid)
                xsec = np.concatenate(xsec)

            return wngrid, xsec

    def compute_xsec_parallel(self,wngrid,T,P,with_progress=True, chunksize=10000,
                     wing_cutoff=25.0, threshold=1e-34, max_workers=4, max_jobs=100):
        import numpy as np
        import concurrent.futures
        import itertools
        from tqdm import tqdm
        import math
        bar = None
        _wngrid = np.sort(wngrid)


        itera = self._linelist.transitions(_wngrid,T, P,wing_cutoff=wing_cutoff,
                                                               chunksize=chunksize,
                                                               threshold=threshold)


        itera = itera

        xsec = []

        out_grid = []
        if self._compute_voigt:
            xsec = np.zeros_like(wngrid)
        
        from tqdm import tqdm
        with tqdm() as t:
            with concurrent.futures.ProcessPoolExecutor() as executor:

                # Schedule the first N futures.  We don't want to schedule them all
                # at once, to avoid consuming excessive amounts of memory.

                futures = {
                    executor.submit(parallel_voigt, task, wing_cutoff, wngrid)
                    for task in itertools.islice(itera, max_jobs)
                }

                while futures:
                    # Wait for the next future to complete.
                    done, futures = concurrent.futures.wait(
                        futures, return_when=concurrent.futures.FIRST_COMPLETED
                    )

                    for fut in done:
                        res,s,e, count = fut.result()
                        
                        if res is None:
                            continue
                        t.update(count)
                        xsec[s:e] += res
                        

                    # Schedule the next set of futures.  We don't want more than N futures
                    # in the pool at a time, to keep memory consumption down.
                    for task in itertools.islice(itera, len(done)):
                        futures.add(
                            executor.submit(parallel_voigt, task, wing_cutoff, wngrid)
                        )

            return wngrid, xsec