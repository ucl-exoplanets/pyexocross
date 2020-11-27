from .exomolstate import ExomolStates

class ExomolTransitionReader:

    @classmethod
    def read_transitions(cls, filename, chunksize=10000):
        import pandas as pd
        df_chunk = pd.read_csv(filename, delim_whitespace=True,usecols=[0,1,2],chunksize=chunksize,header=None)

        for chunk in df_chunk:
            
            chunk.iloc[:,0] -= 1
            chunk.iloc[:,1] -= 1

            yield chunk