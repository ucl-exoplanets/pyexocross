from taurex.log import Logger
class Linelist(Logger):
    
    def __init__(self):
        super().__init__(self.__class__.__name__)
        self._broadeners = {}

    
    def add_broadener(self, broadener, ratio=None):
        species = broadener.species
        if ratio:
            broadener.ratio = ratio
        if species not in self._broadeners:
            self._broadeners[species] = broadener
        self.normalize_broadeners()
    def set_broadener_ratio(self, species, ratio):
        self._broadeners[species].ratio = ratio
        self.normalize_broadeners()
    
    def remove_broadener(self, species):
        if species in self._broadeners:
            del self._broadeners[species]
        


    def normalize_broadeners(self):
        total = sum([v.ratio for v in self._broadeners.values()])
        for v in self._broadeners.values():
            v.ratio/=total

    @property
    def totalTransitions(self):
        raise NotImplementedError

    @property
    def molecularMass(self):
        raise NotImplementedError


    def compute_doppler(self, temperature, df):
        from .util import doppler_broad
        import math
        freq = df['v_if'].values
        return doppler_broad(freq, self.molecularMass, temperature)

    def compute_intensity(self,trans, temperature,pf=None, wn_filter=None):
        from .constants import PLANCK, PI, KBOLTZ, SPDLIGT, c2
        import numexpr as ne

        if wn_filter is None:
            wn_filter = slice(None)

        gtot_f = trans["g_tot'"].values[wn_filter]
        Aif_v = trans['A_if'].values[wn_filter]
        E_i = trans['E"'].values[wn_filter]
        v = trans['v_if'].values[wn_filter]
        if pf is None:
            pf = self.compute_partition(temperature, trans)
        elif not isinstance(pf, float):
            pf = pf.Q(temperature)

        T = temperature

        return ne.evaluate('gtot_f*Aif_v*exp(-c2*E_i/T)*(1-exp(-c2*v/T))/(8*PI*SPDLIGT*v*v*pf)')
    
    def get_transitions(self,min_wn, max_wn, chunksize=10000):
        raise NotImplementedError

    
    def estimated_count(self, wngrid):
        return self.totalTransitions

    def transitions(self, wngrid, temperature, pressure,pf=None, wing_cutoff=25.0, chunksize=10000, threshold=1e-34):
        import numpy as np
        self.normalize_broadeners()
        for v in self._broadeners.values():
            v.set_temperature_pressure(temperature, pressure)
        
        min_wn, max_wn = wngrid.min()-wing_cutoff, wngrid.max()+wing_cutoff


        for df in self.get_transitions(min_wn, max_wn, chunksize=chunksize):
            read_chunks = len(df)
            v = df['v_if'].values
            transition_filter = (v >= min_wn) & (v <= max_wn)

            df = df[transition_filter]
            if len(df) == 0:
                yield None, None, None, None, read_chunks  
                continue      
            I = self.compute_intensity(df, temperature, pf=pf)
            threshold_filter = I >= threshold
            I = I[threshold_filter]
            df = df[threshold_filter]
            if len(self._broadeners) > 0:
                gamma = sum([g.gamma(df) for g in self._broadeners.values()])
            else:
                gamma = np.zeros_like(I)
            v = df['v_if'].values
            I = I*self.get_abundance(df)
            gamma = gamma
            doppler = self.compute_doppler(temperature,df)
            yield v, I, gamma, doppler, read_chunks

    def get_abundance(self, df):
        return 1.0


