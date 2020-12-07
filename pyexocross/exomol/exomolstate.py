import pandas as pd
from ..constants import PLANCK, PI, KBOLTZ, SPDLIGT, c2

class ExomolStates:

    def __init__(self,filename, fwf_definition):

        self._df = self._read_state(filename, *fwf_definition)
        self.find_quanta_start()

    def find_quanta_start(self):
        self._start = 4
        if 'lftime' in self._df.columns:
            self._start += 1
        if 'lande-g' in self._df.columns:
            self._start += 1

    @property
    def pandasFrame(self):
        return self._df

    def _read_state(self, filename, header, width, dtype):
        import pandas as pd
        import numpy as np
        df = pd.read_csv(filename,delim_whitespace=True,names=header)
        
        for key, value in dtype.items():
            if value in (np.int64, np.float64, ):
                d = value

                df[key] = pd.to_numeric(df[key], errors='coerce')
                if value is np.int64:
                    d = 'Int64'
                    df[key] = df[key].astype(d)
        return df

    def Q(self, T):
        import numexpr as ne
        E = self._df.E.values
        gns = self._df.g_tot.values

        return ne.evaluate('sum(gns*exp(-c2*E/T))')


    def transition_states(self, transitions):
        import numexpr as ne
        idf = transitions.iloc[:,0]
        idi = transitions.iloc[:,1]
        upper = self._df.iloc[idf]
        lower = self._df.iloc[idi]
        # if pf is None:
        #     pf = self.Q(temperature)
        # elif not isinstance(pf, float):
        #     pf = pf.Q(temperature)
            


        Aif = transitions.iloc[:,2]

        transition_frame = {}

        v = upper.E.values - lower.E.values
        Aif_v = Aif.values
        gtot_f = upper.g_tot.values
        E_i = lower.E.values

        transition_frame['v_if'] = v
        
        # transition_frame['Iif'] = \
        #     ne.evaluate('gtot_f*Aif_v*exp(-c2*E_i/T)*(1-exp(-c2*v/T))/(8*PI*SPDLIGT*v*v*pf)')

        transition_frame['A_if'] = Aif_v
        



        transition_frame["E'"] = upper.E.values
        transition_frame["g_tot'"] = gtot_f
        transition_frame["J'"] = upper.J.values
        for col in upper.columns[self._start:]:
            transition_frame[f"{col}'"] = getattr(upper,col).values

        transition_frame['E"'] = E_i
        transition_frame['g_tot"'] = lower.g_tot.values
        transition_frame['J"'] = lower.J.values
        for col in lower.columns[self._start:]:
            transition_frame[f'{col}"'] = getattr(lower,col).values

        df=pd.DataFrame(transition_frame)

        # if threshold is not None:
        #     df = df[df.Iif >= threshold]
        
        return df









