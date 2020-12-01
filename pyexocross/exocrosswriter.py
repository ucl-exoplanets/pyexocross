from taurex.log import Logger
from taurex.util.util import calculate_weight, sanitize_molecule_string, \
    conversion_factor
from .util import convert_to_wavenumber
import numpy as np
import os
import glob
import tabulate

class BroadenerInput:

    def __init__(self, molecule, gamma, n, t0=298.0, p0=1.0, ratio=0.5,
                 filename=None, broadener_type=None, delta=None):
        self.molecule = molecule
        self.gamma = gamma
        self.n = n
        self.t0 = t0
        self.p0 = p0
        self.ratio = ratio
        self.filename = filename
        self.broadener_type = broadener_type
        self.delta = None

    def generate_input(self, path='.'):
        output_string = []
        output_string.append(f'{self.molecule} gamma {self.gamma}  n {self.n} t0 {self.t0}')
        if self.filename is not None:
            output_string.append(f'file {os.path.join(path,self.filename)}')
            if self.broadener_type is not None:
                output_string.append(f'model {self.broadener_type}')
            else:
                raise ValueError(f'No model specified for broadener {self.molecule}')

        if self.delta is not None:
            output_string.append(f'delta {self.delta}')

        output_string.append(f'ratio {self.ratio}')

        return ' '.join(output_string)


class ExocrossInput(Logger):

    def __init__(self, molecule=None, temperature=500.0,
                 pressure=1.0, linelist=None, Q=None, state_file=None,
                 file_prefix=None, path='.'):
        super().__init__(self.__class__.__name__)
        self.temperature = 500.0
        self.pressure = conversion_factor('Pa', 'bar') * pressure
        self.molecule = molecule
        self.linelist = None
        self.partition_function = Q
        self._broadeners = {}
        self.determine_mass()
        self.offset = 25.0
        self._npoints = 1001
        self.threshold = 0.0
        self._range = (0.001, 10000.0)
        self._state_file = state_file
        self._transitions = []
        self._prefix = None
        self.set_file_path(path)
        
        self.set_file_prefix(file_prefix)

    
    def set_file_path(self, path, reloadp=True):
        self._path = path
        if self._prefix is not None and reloadp:
            self.discover_all()


    def set_file_prefix(self, prefix, reloadp=True):
        self._prefix = prefix
        if self._path is not None and reloadp:
            self.discover_all()
    
    def discover_all(self):
        self.discover_state_file()
        self.discover_transitions()
        self.discover_partition()

    def discover_state_file(self):
        if self._prefix is not None and self._path is not None:
            state_file = f'{self._prefix}.states'
            filename = os.path.join(self._path, state_file)
            if os.path.isfile(filename):
                self.info(f'Discovered state_file {state_file}')
                self.set_state_file(filename)

    def discover_transitions(self):
        import pathlib
        if self._prefix is not None and \
           self._path is not None:
            trans_filename = f'{self._prefix}__*.trans'
            filename = os.path.join(self._path, trans_filename)

            trans_list = glob.glob(filename)
            just_file = [pathlib.Path(f).stem for f in trans_list]
            wavenumbers = [f[-11:].split('-') for f in just_file]
            float_wn = [(float(mn), float(mx)) for mn, mx in wavenumbers]

            self._transitions = list(zip(float_wn, trans_list))

            self._transitions.sort(key=lambda x: x[0])

            self.info('Discovered Transitions \n\n %s \n\n',
                      tabulate.tabulate(self._transitions,
                                        ['Spectral Range', 'Filename']))


    def discover_partition(self):
        if self._prefix is not None and self._path is not None:
            state_file = f'{self._prefix}.pf'
            filename = os.path.join(self._path, state_file)
            if os.path.isfile(filename):
                self.info(f'Discovered partition file {state_file}')
                self.load_pf(filename)

    def load_pf(self, filename):
        from scipy.interpolate import interp1d
        self._pf_func = None
        if filename is not None:
            with open(filename, 'r') as f:
                pf_array = np.loadtxt(f)
            self._pf_func = interp1d(pf_array[:, 0], pf_array[:, 1])

    @property
    def Q(self):
        if self.partition_function is not None:
            return self.partition_function
        elif self._pf_func is not None:
            return self._pf_func(self.temperature)[()]
        else:
            return None

    def determine_mass(self):
        self.mass = calculate_weight(sanitize_molecule_string(self.molecule))
        self.info('Auto determined mass as %.3f', self.mass)

    def normalize_broadeners(self):

        total_broad = sum([b.ratio for b in self._broadeners.values()])

        self.info('Total broadeners sum to %.3f', total_broad)

        for b in self._broadeners.values():
            b.ratio /= total_broad

    def add_broadener(self, broadener):
        molecule = broadener.molecule
        if molecule not in self._broadeners:
            self._broadeners[molecule] = broadener
        else:
            raise KeyError(f'Broadener {molecule} already exists')

    def remove_broadener(self, molecule):
        if molecule not in self._broadeners:
            raise KeyError(f'Broadener {molecule} does not exist')
        else:
            del self._broadeners[molecule]

    def get_broadener(self, molecule):
        if molecule not in self._broadeners:
            raise KeyError(f'Broadener {molecule} does not exist')
        else:
            return self._broadeners[molecule]

    def set_broadener_ratio(self, molecule, ratio):
        broad = self.get_broadener(molecule)
        broad.ratio = ratio

    def set_state_file(self, state_file):
        self._state_file = state_file

    def set_pressure(self, pressure, unit='Pa'):
        factor = conversion_factor(unit, 'bar')
        self.pressure = pressure*factor

    def set_molar_mass(self, mass):
        self.mass = mass

    def set_range(self, spectral_range, units='k'):

        wm_range = convert_to_wavenumber(spectral_range, units)

        self.range = (min(wm_range), max(wm_range))

    def valid_transitions(self):

        start, end = self.range

        return [filename for (min_wn, max_wn), filename in self._transitions 
                if min_wn <= end and max_wn >= start]

    @property
    def Npoints(self):
        return self._npoints

    @Npoints.setter
    def Npoints(self, value):
        int_value = int(abs(value))
        # Determine oddness
        if int_value & 1 == 0:
            int_value += 1
        self._npoints = int_value

    def add_transition(self, transition_file):
        self._transitions.append(transition_file)

    def add_transition_list(self, transition_list):

        self._transitions.extend(transition_list)

    def clear_transitions(self):
        self._transitions = []

    def generate_input(self, output_filename):

        output_string = []

        output_string.append(f'Temperature {self.temperature} (K)')
        output_string.append(f'Range {self.range[0]} {self.range[1]} (cm-1)')
        output_string.append(f'Npoints {self.Npoints}')

        output_string.append(f'pressure {self.pressure} (bar)')

        output_string.append('absorption')
        output_string.append('voigt')

        output_string.append(f'threshold {self.threshold}')

        output_string.append(f'mass {self.mass}')
        Q = self.Q
        if Q is not None:
            output_string.append(f'PF {Q}')

        if len(self._broadeners) > 0:
            self.normalize_broadeners()
            output_string.append('\nspecies')
            for b in self._broadeners.values():
                output_string.append(f'  {b.generate_input(path=self._path)}')
            output_string.append('end\n')

        output_string.append(f'\nStates {self._state_file}')
        output_string.append(f'\nTransitions')
        for t in self.valid_transitions():
            output_string.append(f'  {t}')
        output_string.append('end\n')

        output_string.append(f'Output {output_filename}')

        return '\n'.join(output_string)
