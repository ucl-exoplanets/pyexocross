from .linelist import Linelist
from .voigt_functions import Voigt



def parallel_voigt(args, wing_cutoff=25.0, wngrid=None, accurate=False, function='scipy'):
    v, I, gamma, doppler,count = args

    if v is None or len(v)==0:
        return None,None,None,count
    voigt = Voigt(voigt_function=function, accurate_sum=accurate)
    result = voigt.voigt(wngrid, v, I, doppler, gamma,cutoff=wing_cutoff)
    return result[0],result[1],result[2],count 

def create_jobs(linelist_iterator, wing_cutoff, wngrid, queue, num_workers, accurate=False, function='scipy'):
    import concurrent.futures
    import numpy as np
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        #for task in linelist_iterator
        for v, I, gamma, doppler,count in linelist_iterator:
            if v is None or len(v) == 0:
                queue.put(executor.submit(parallel_voigt, (v, I, gamma, doppler,count), wing_cutoff, wngrid,
                                                            accurate, function))
                continue
            v_s = np.array_split(v,num_workers)
            I_s = np.array_split(I,num_workers)
            gamma_s = np.array_split(gamma,num_workers)
            doppler_s = np.array_split(doppler,num_workers)
            count_s = [count//num_workers for x in v_s]
            for task in zip(v_s,I_s, gamma_s, doppler_s, count_s):
                # if count_s == 0:
                #     continue
                queue.put(executor.submit(parallel_voigt, task, wing_cutoff, wngrid, accurate, function))
    queue.put(False)

class PyExocross:

    def __init__(self, linelist: Linelist, compute_voigt=True):
        self._linelist = linelist
        self._voigt = Voigt()
    
        self._compute_voigt = compute_voigt
        self.set_voigt_function('scipy')
        self.use_accurate_sum(False)

    def set_voigt_function(self, voigt_function):
        self._voigt_func = voigt_function
        self._voigt.set_voigt_function(voigt_function)

    def use_accurate_sum(self, accurate_sum):
        self._voigt_sum = accurate_sum
        self._voigt.use_accurate_sum(accurate_sum)

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
        total_size = self._linelist.estimated_count(wngrid)
        with tqdm(total=total_size) as t:
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
                     wing_cutoff=25.0, threshold=1e-34, max_workers=4, max_jobs=None):
        import numpy as np
        
        import itertools
        from tqdm import tqdm
        import math
        import threading
        import queue
        bar = None
        _wngrid = np.sort(wngrid)
        if max_jobs is None:
            max_jobs = max_workers*4

        itera = self._linelist.transitions(_wngrid,T, P,wing_cutoff=wing_cutoff,
                                                               chunksize=chunksize,
                                                               threshold=threshold)


        xsec = []

        out_grid = []
        if self._compute_voigt:
            xsec = np.zeros_like(wngrid)
        
        from tqdm import tqdm
        total_size = self._linelist.estimated_count(wngrid)
        job_queue = queue.Queue(max_jobs)
        job_creator = threading.Thread(target=create_jobs,
                                        args=(itera,
                                                wing_cutoff, wngrid, job_queue, max_workers,self._voigt_sum,self._voigt_func))
        job_creator.start()
        with tqdm(total=total_size) as t:
            while True:
                future = job_queue.get()
                if not future:
                    job_queue.task_done()
                    break
                job_queue.task_done()
                res, s, e, count = future.result()
                if count is not None:
                    t.update(count)
                if res is None:
                    continue
                
                xsec[s:e]+=res
        job_creator.join()

        return wngrid, xsec