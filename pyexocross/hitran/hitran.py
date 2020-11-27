import numpy as np
from ..linelist import Linelist
from ..broadener import Broadener
import os


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

    def __init__(self,filename,
                 ):
        super().__init__()
        self._filename = filename
        self.load_hitran_file(filename)
        filesize = os.path.getsize(filename)
        self._total_transitions = filesize/160

    @property
    def totalTransitions(self):
        return self._total_transitions

    def add_self_broadener(self):
        if 'self' not in self._broadeners:
            self.add_broadener(HITRANSelfBroadener())

    def add_air_broadener(self):
        if 'air' not in self._broadeners:
            self.add_broadener(HITRANAirBroadener())

    def load_hitran_file(self, filename):
        with open(filename,'r') as f:
            line = f.readline()
            self._molid = int(line[:2])
            self._isoid = int(line[2:3])

    def compute_partition(self, temperature):
        from .hapi import partitionSum
        return partitionSum(self._molid, self._isoid, temperature)



    @property
    def molecularMass(self):
        from .hapi import molecularMass
        return molecularMass(self._molid, self._isoid)


    def get_transitions(self, chunksize=10000):
        from .util import read_hitran_pandas
        yield from read_hitran_pandas(self._filename, chunksize=chunksize)







    