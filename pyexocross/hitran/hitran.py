import numpy as np
from ..linelist import Linelist
class HITRANLinelist(Linelist):

    def __init__(self,filename,
                 ):
        self._filename = filename
        # self.load_hitran_file(filename)
        self.set_broadener_ratio(0.0,1.0)
    def set_broadener_ratio(self, self_ratio, air_ratio):
        self._self_ratio = self_ratio
        self._air_ratio = air_ratio

        total = self._self_ratio + self._air_ratio
        self._self_ratio /= total
        self._air_ratio /= total
    




    def hitran_iterator(self, filename):
        from .util import read_hitran
        self._transition_df
        # molid, isoid, v, S, A, gamma_air, gamma_self, E_lower, n_air, \
        # delta_air, upper_quant, lower_quant,gf,gi = \
        #     zip(*[(mol, isom, vif, Sif, Aif, gair, gself, El, nair, delta_air, upq, lq,gf,gi) for \
        #           mol, isom, vif, Sif, Aif, gair, gself, El, nair, delta_air, upq, lq,gf,gi in read_hitran(filename)])

        # self._molid = molid[0]
        # self._isoid = isoid[0]
        # self._v = np.array(v)
        # self._A = np.array(A)
        # self._gamma_air = np.array(gamma_air)
        # self._gamma_self = np.array(gamma_self)
        # self._n_air = np.array(n_air)
        # self._E_lower = np.array(E_lower)
        # self._delta_air = np.array(delta_air)
        # self._gf = np.array(gf)
        # self._gi = np.array(gi)




    @property
    def totalTransitions(self):
        return len(self._v)

    def compute_partition(self, temperature):
        from .hapi import partitionSum
        return partitionSum(self._molid, self._isoid, temperature)
    
    def compute_intensity(self, temperature,wn_filter=None):
        from ..constants import PLANCK, PI, KBOLTZ, SPDLIGT, c2
        import numexpr as ne

        if wn_filter is None:
            wn_filter = slice(None)

        gtot_f = self._gf[wn_filter]
        Aif_v = self._A[wn_filter]
        E_i = self._E_lower[wn_filter]
        v = self._v[wn_filter]
        pf = self.compute_partition(temperature)
        T = temperature

        return ne.evaluate('gtot_f*Aif_v*exp(-c2*E_i/T)*(1-exp(-c2*v/T))/(8*PI*SPDLIGT*v*v*pf)')

    def compute_gamma(self, temperature, pressure, transition_filter=None):
        t0 = 296.0
        p0 = 1.0
        T = temperature
        P = pressure
        gamma_air = 0.0
        if transition_filter is None:
            transition_filter = slice(None)

        if self._air_ratio > 0.0:
            gamma_air = self._gamma_air[transition_filter]*((t0/T)**self._n_air[transition_filter])*(P/p0)*self._air_ratio
        
        gamma_self = 0.0
        # if self._self_ratio > 0.0:
        #     gamma_self = self._gamma_self[transition_filter]*((t0/T)**self._n_self[transition_filter])*(P/p0)*self._self_ratio

        return gamma_air + gamma_self

    @property
    def molecularMass(self):
        from .hapi import molecularMass
        return molecularMass(self._molid, self._isoid)


    def compute_doppler(self, temperature, freq):
        from ..constants import KBOLTZ, SPDLIGT
        import math
        return math.sqrt(2*KBOLTZ*math.log(2)/self.molecularMass)*freq/SPDLIGT


    
    def transitions(self, wngrid, temperature, pressure, wing_cutoff=25.0, chunksize=10000, threshold=1e-34):
        
        min_wn, max_wn = wngrid.min()-wing_cutoff, wngrid.max()+wing_cutoff

        transition_filter = (self._v >= min_wn) & (self._v <= max_wn)
        v = self._v[transition_filter]
        I = self.compute_intensity(temperature, wn_filter=transition_filter)
        gamma = self.compute_gamma(temperature, pressure, transition_filter)

        threshold_filter = I >= threshold

        v = v[threshold_filter]
        I = I[threshold_filter]
        gamma = gamma[threshold_filter]
        doppler = self.compute_doppler(temperature,v)

        total_transitions = len(v)

        for x in range(0, len(v), chunksize):
            yield v[x:x+chunksize], I[x:x+chunksize], gamma[x:x+chunksize], doppler[x:x+chunksize]







    