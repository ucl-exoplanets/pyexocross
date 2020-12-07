from ..linelist import Linelist
import glob
import os
import tabulate


class ExomolLinelist(Linelist):

    def __init__(self, path_to_linelist=None):
        super().__init__()
        self._prefix = None
        self._path = None
        self.set_file_path(path_to_linelist, reloadp=True)
        #if path_to_linelist is not None:
            #def_file


    def set_file_path(self, path, reloadp=True):
        self._path = path
        if reloadp:
            self.discover_def()

    @property
    def totalTransitions(self):
        return self._exodef._total_transitions

    def set_file_prefix(self, prefix, reloadp=True):
        self._prefix = prefix
        if self._path is not None and reloadp:
            self.discover_all()
    
    def discover_all(self):
        self.discover_state_file()
        self.discover_transitions()
        self.discover_partition()
        self.discover_broadeners()

    def discover_def(self):
        
        from .exomoldef import ExomolDef
        if self._path is not None:
            def_file = glob.glob(os.path.join(self._path,'*.def'))[0]
            self._exodef = ExomolDef(def_file)
            self.set_file_prefix(self._exodef.filePrefix,reloadp=True)

    def discover_state_file(self):
        if self._prefix is not None and self._path is not None:
            state_file = f'{self._prefix}.states'
            filename = os.path.join(self._path, state_file)
            if os.path.isfile(filename):
                self.info(f'Discovered state_file {state_file}')
                self.set_state_file(filename)
            else:
                state_file = f'{self._prefix}.states.bz2'
                filename = os.path.join(self._path, state_file)
                if os.path.isfile(filename):
                    self.info(f'Discovered state_file {state_file}')
                    self.set_state_file(filename)


    def discover_transitions(self):
        import pathlib
        import numpy as np
        if self._prefix is not None and \
           self._path is not None:
            trans_filename = f'{self._prefix}*.trans'
            trans_filename_bz2 = f'{self._prefix}*.trans.bz2'
            filename = os.path.join(self._path, trans_filename)
            trans_list = glob.glob(filename)
            if len(trans_list) == 0:
                filename = os.path.join(self._path, trans_filename_bz2)
                trans_list = glob.glob(filename)
            if len(trans_list) == 1:
                filename = os.path.join(self._path, trans_list[0])
                self._transitions = [((0.0, np.inf),filename)]
            else:      

                just_file = [pathlib.Path(f).stem for f in trans_list]
                just_file = [f if not f.endswith('.trans') else f[:-6] for f in just_file]

                wavenumbers = [f[-11:].split('-') for f in just_file]


                float_wn = [(float(mn), float(mx)) for mn, mx in wavenumbers]

                self._transitions = list(zip(float_wn, trans_list))

                self._transitions.sort(key=lambda x: x[0])

            self.info('Discovered Transitions \n\n %s \n\n',
                      tabulate.tabulate(self._transitions,
                                        ['Spectral Range', 'Filename']))


    def compute_partition(self, temperature, df):
        return self._pf.Q(temperature)

    def discover_partition(self):
        from .exomolpf import ExomolPF
        if self._prefix is not None and self._path is not None:
            state_file = f'{self._prefix}.pf'
            filename = os.path.join(self._path, state_file)
            if os.path.isfile(filename):
                self.info(f'Discovered partition file {state_file}')
                self._pf = ExomolPF(filename)

    @property
    def molecularMass(self):
        return self._exodef._mass


    def discover_broadeners(self):
        from pathlib import Path
        self._avail_broadeners = []
        if self._prefix is not None and \
           self._path is not None:
            broad_filename = f'{self._exodef._molecule_slug}__*.broad'
            filename = os.path.join(self._path, broad_filename)

            broad_list = glob.glob(filename)

            files = {Path(b).stem.split('_')[-1]: b for b in broad_list}

            self._avail_broadeners = files

    @property
    def availableBroadeners(self):
        return tuple(self._avail_broadeners.keys())


    def add_available_broadener(self, broad_name, ratio=1.0):
        if broad_name in self.availableBroadeners:
            filename = self._avail_broadeners[broad_name]
            if broad_name in self._exodef.availableBroadeners:
                broad = self._exodef._broadener_defs[broad_name].generate_exomolbroadener(filename)
                self.add_broadener(broad,ratio=ratio)
            else:
                from .exomolbroads import ExomolBroadener
                gamma, n =self._exodef._default_gamma, self._exodef._default_n
                broad = ExomolBroadener(gamma, n, filename=filename)
                self.add_broadener(broad,ratio=ratio)
    

    def add_default_broadener(self, ratio=1.0):
        from .exomolbroads import ExomolBroadener
        if 'default' not in self.availableBroadeners:
            gamma, n =self._exodef._default_gamma, self._exodef._default_n
            broad = ExomolBroadener(gamma, n, species='default')
            self.add_broadener(broad,ratio=ratio)
            

    def valid_transitions(self, start, end):

        return [filename for (min_wn, max_wn), filename in self._transitions 
                if min_wn <= end and max_wn >= start]

    def set_state_file(self, state_file):
        
        self._state_file = state_file
        self._state = self._exodef.read_state(state_file)


    def get_transitions(self,min_wn, max_wn, chunksize=10000):
        from .exomoltransition import ExomolTransitionReader as etr
        files = self.valid_transitions(min_wn, max_wn)

        for f in files:
            for c in etr.read_transitions(f,chunksize=chunksize):
                yield self._state.transition_states(c)
