import numpy as np
from ..linelist import Linelist
from ..broadener import Broadener
import pathlib


class HITRANSelfBroadener(Broadener):

    def __init__(self, ratio=1.0):
        super().__init__(ratio=ratio)
    
    @property
    def species(self):
        return 'self'

    def calculate_gamma(self, transitions):
        from ..util import compute_gamma
        return compute_gamma(transitions['gamma_self'].values,1.0, self.T, self.P,
                             296.0, 1.0)

class HITRANAirBroadener(Broadener):

    def __init__(self, ratio=1.0):
        super().__init__(ratio=ratio)
    
    @property
    def species(self):
        return 'air'

    def calculate_gamma(self, transitions):
        from ..util import compute_gamma
        return compute_gamma(transitions['gamma_air'].values,transitions['n_air'].values, self.T, self.P,
                             296.0, 1.0)


class HITRANLinelist(Linelist):

    def __init__(self,filename, iso_abundances=None
                 ):
        super().__init__()
        self._filename = filename
        self.load_hitran_file(filename)
        filesize = pathlib.Path(filename).stat().st_size
        self._total_transitions = filesize//(self._total_line)
        
    @property
    def totalTransitions(self):
        return self._total_transitions

    def add_self_broadener(self, ratio=1.0):
        if 'self' not in self._broadeners:
            self.add_broadener(HITRANSelfBroadener(ratio=ratio))

    def add_air_broadener(self, ratio=1.0):
        if 'air' not in self._broadeners:
            self.add_broadener(HITRANAirBroadener(ratio=ratio))

    def load_hitran_file(self, filename):
        from .hapi import molecularMass
        with open(filename,'r') as f:
            line = f.readline()
            self._molid = int(line[:2])
            self._total_line = len(line)+1
        self.discover_iso()


    def discover_iso(self):
        from .hapi import ISO, molecularMass, abundance
        self._isotopalogues = np.array([k[1] for k in ISO if k[0] == self._molid],dtype=np.int)
        max_iso = self._isotopalogues.max()
        self._molmasses = np.empty(shape=(max_iso))
        self._abundance_vals = np.empty(shape=(max_iso))
        for iso in self._isotopalogues:
            self._molmasses[iso-1] = molecularMass(self._molid,iso)
            self._abundance_vals[iso-1] = abundance(self._molid, iso)
        self._abundance_vals/=self._abundance_vals.sum()

    @property
    def molecule(self):
        from .hapi import moleculeName
        return moleculeName(self._molid)

        
    def compute_partition(self, temperature, df):
        from .hapi import partitionSum
        
        isoid = df['IsoID'].values -1
        return self._iso_partition[isoid]


    def compute_doppler(self, temperature, df):
        from ..util import doppler_broad
        import math
        freq = df['v_if'].values
        masses = self._molmasses[df['IsoID'].values-1]
        return doppler_broad(freq, masses, temperature)

    def get_transitions(self,min_wn, max_wn, chunksize=10000):
        from .util import read_hitran_pandas
        yield from read_hitran_pandas(self._filename, chunksize=chunksize)

    def get_abundance(self, df):
        abundances = self._abundance_vals[df['IsoID'].values-1]
        return abundances


    def transitions(self, wngrid, temperature, pressure,pf=None, wing_cutoff=25.0, chunksize=10000, threshold=1e-34):
        from .hapi import partitionSum
        self._iso_partition = np.array([partitionSum(self._molid, x+1, temperature) 
                                        for x in range(self._molmasses.shape[0])])
        yield from super().transitions(wngrid, temperature, pressure,pf, wing_cutoff, chunksize, threshold)