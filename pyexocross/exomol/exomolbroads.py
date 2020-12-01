import os
import pathlib
import numpy as np
from ..broadener import Broadener
def _filter_func(trans, label, value):
    return trans[label] == value



class ExomolBroadener(Broadener):

    def __init__(self, default_gamma0, default_n0,t0=298.0, p0=1.0,
                 label_defs={},filename=None, ratio=1.0):
        super().__init__(ratio=ratio)
        
        self._t0 = t0
        self._p0 = p0

        self._ratio = 1.0

        self._gamma0 = default_gamma0
        self._n0 = default_n0

        

        self._label_functions = {}
        self._broadener_functions = []
        self._broadener_values = []
        self._precomputed_values = []
        self._default_gamma = 0.0
        for k,v in label_defs.items():
            label = k
            self.add_new_broad_label(label,v)
        
        if filename is not None:
            self.load_file(filename)
    
    def load_file(self, filename):
        if not os.path.isfile(filename):
            raise FileNotFoundError(f'No file with path {filename} exists')
        self._species = pathlib.Path(filename).stem.split('_')[-1].strip()
        with open(filename, 'r') as f:

            for line in f:
                label, gamma, n, *quanta = line.split()
                gamma = float(gamma)
                n = float(n)
                quanta = [int(q) for q in quanta]
                if self.label_exists(label):
                    func = self._label_functions[label]
                    self._broadener_functions.append(lambda trans, args=quanta, func=func: func(trans,*args))
                    self._broadener_values.append((gamma,n))

    def set_temperature_pressure(self, T, P):
        super().set_temperature_pressure(T, P)
        self._precomputed_values = [self.compute_gamma(gamma,n,T,P) for gamma, n in self._broadener_values]
        self._default_gamma = self.compute_default_gamma(T,P)


    def label_exists(self, label):
        return label in self._label_functions

    def add_new_broad_label(self, label, label_definition):
        from operator import and_
        from functools import reduce
        if label == 'a0':
            self._label_functions['a0'] = lambda trans, *args: _filter_func(trans,'J"',args[0])
        else:
            self._label_functions[label] = lambda trans,*args,label_def=label_definition: np.logical_and(*[_filter_func(trans,l,a) for l,a in zip(label_def,args)])
            #self._label_functions[label] = lambda trans, **args: all(_)
    
    def compute_gamma(self, gamma0,n0, T, P):
        
        return gamma0*((self._t0/T)**n0)*(P/self._p0)*self.ratio

    def compute_default_gamma(self, T,P):
        return self.compute_gamma(self._gamma0, self._n0, T, P)

    @property
    def species(self):
        return self._species

    @property
    def ratio(self):
        return self._ratio
    
    @ratio.setter
    def ratio(self, value):
        self._ratio = value

    
    def gamma(self, transitions):
        import numpy as np

        gamma = np.zeros(shape=(len(transitions)))
        filt = np.empty(shape=(len(transitions)),dtype=np.bool)
        filt[...] = True
        for broad,gamma_val in zip(self._broadener_functions,self._precomputed_values):
            # Get applicable gamma
            broad_filt = broad(transitions)
            # See leftover transitions
            applic = filt & broad_filt
            # If we have anny then apply it
            if np.any(applic):
                gamma[applic] = gamma_val
            # Now exclude them from the search
            filt &= ~broad_filt
        # If we have any leftovers, apply the default
        if np.any(filt):
            gamma[filt] = self._default_gamma

        return gamma
            


