from taurex.log import Logger


class LinesReader:

    def __init__(self, lines):
        self._lines = lines
        self._count = 0
    
    def skip(self, num=1):
        self._count += num

    def read_int(self, skip=1):
        val = int(self._lines[self._count])
        self.skip(skip)
        return val
    
    def read_float(self, skip=1):
        val = float(self._lines[self._count])
        self.skip(skip)
        return val
    
    def read_float_array(self, skip=1):
        line = self.read_string()
        split = line.split()
        return [float(s) for s in split]
    
    def read_string(self, skip=1):
        val = self._lines[self._count]
        self.skip(skip)
        return val
    
    def read_bool(self, skip=1):
        val = int(self._lines[self._count])
        
        self.skip(skip)
        return val == 1

    def reset():
        self._count = 0

class BroadenerData:

    def __init__(self, molecule, filename, Jmax, default_gamma, default_n):
        self._molecule = molecule.strip()
        self._filename = filename
        self._default_gamma = default_gamma
        self._default_n = default_n
        self._Jmax = Jmax
        self._avail_codes = []
        self._quanta={}
    def add_code(self, quanta_code, quanta):
        self._avail_codes.append(quanta_code)
        quanta.insert(0,'J"')
        self._quanta[quanta_code] = quanta

    @property
    def molecule(self):
        return self._molecule

    @property
    def availableCodes(self):
        return self._avail_codes

    def generate_input(self, maximum_model='JJ', broadener_path='.'):
        from .exocrosswriter import BroadenerInput
        import os
        bb = BroadenerInput(self._molecule, self._default_gamma,
                            self._default_n, filename=os.path.join(broadener_path, self._filename),
                            broadener_type='JJ' if 'a1' in self._avail_codes else 'J')
        return bb
    
    def generate_exomolbroadener(self, filename=None):
        from .exomolbroads import ExomolBroadener
        return ExomolBroadener(self._default_gamma, self._default_n,label_defs=self._quanta,filename=filename)

class ExomolDef(Logger):

    def __init__(self, exomol_def_file):
        super().__init__(self.__class__.__name__)

        self.info(f'Opening {exomol_def_file}')
        with open(exomol_def_file, 'r') as f:
            unclean_exocross_lines = f.read().splitlines()
            self.exocross_lines = [s.split('#')[0].strip()
                                   for s in unclean_exocross_lines]

        self.parse_definition()

    def parse_definition(self):
        
        lr = LinesReader(self.exocross_lines)

        if lr.read_string() != 'EXOMOL.def':
            raise IOError('Incorrect EXOMOL def header')

        lr.skip(1)
        
        self._molecule_slug = lr.read_string()

        self._linelist_name = lr.read_string()

        self._version_number = lr.read_string()

        self._inchikey = lr.read_string()

        self._natoms = lr.read_int()

        self.info(f'Molecule is {self._molecule_slug}')
        self.info(f'Linelist: {self._linelist_name} '
                  f'Version: {self._version_number}')
        while(True):
            try:
                arr = lr.read_float_array()
                self._mass = arr[0]
                test = arr[1]
                break
            except (IndexError, ValueError, ):
                continue
        self._symmetry_group = lr.read_string()

        num_irr = lr.read_int()
        lr.skip(num_irr * 3)
        self._max_temp = lr.read_float()
        self._num_broadeners = lr.read_int()
        self._dipole_avail = lr.read_bool()
        self._no_cross = lr.read_int()
        self._no_ktab = lr.read_int()

        self._life_avail = lr.read_bool()
        self._landeg_avail = lr.read_bool()


        self._num_states = lr.read_int()
        num_cases = lr.read_int()
        self._quanta_cases = {}
        for case in range(num_cases):
            case_label = lr.read_string()
            no_quanta = lr.read_int()
            quanta_definition = []

            for q in range(no_quanta):
                label = lr.read_string()
                
                form =lr.read_string() 
                form = form.split()[1].strip()
                descrp = lr.read_string()

                quanta_definition.append((label, form, descrp))
        

            self._quanta_cases[case_label] = quanta_definition
        
        self._total_transitions = lr.read_int()
        self._num_trans_files = lr.read_int()
        self._max_wavenumber = lr.read_float()
        self._highest_complete = lr.read_float()
        self._max_temp_q = lr.read_float()
        self._t_step = lr.read_float(skip=2)
        self._default_gamma = lr.read_float()
        self._default_n = lr.read_float()
        self._broadener_defs = {}

        if self._num_broadeners > 0:
            for b in range(self._num_broadeners):

                broadener_label = lr.read_string()
                self.info(f'Reading broadener {broadener_label}')
                broadener_filename = lr.read_string()
                jmax = lr.read_int()
                default_gamma = lr.read_float()
                default_n = lr.read_float()

                new_broad = BroadenerData(broadener_label, broadener_filename,
                                          jmax, default_gamma, default_n)
                self._broadener_defs[broadener_label] = new_broad

                n_broad_quanta_set = lr.read_int()
                for x in range(n_broad_quanta_set):
                    code_label = lr.read_string(skip=2)
                    
                    no_quanta = lr.read_int()
                    quanta =[lr.read_string() for x in range(no_quanta)]
                    new_broad.add_code(code_label, quanta)

    def _pandas_state_fwf(self):
        import re
        import numpy as np

        widths = [12,12,6,7]
        headers = ['i', 'E', 'g_tot', 'J']


        if self._life_avail:
            widths.append(12)
            headers.append('lftime')
        
        if self._landeg_avail:
            widths.append(12)
            headers.append('lande-g')
        
        # Let pandas auto determine above types

        # We will be specific about the quanta
        dtype = {}
        form_conv = {'d' : np.int64,
                     'f' : np.float64,
                     's' : str,
                     'i' : np.int64 }
        for case in self._quanta_cases.values():
            for label, form, desc in case:
                headers.append(label)
                wid = re.findall(r'\d+',form)[0]
                widths.append(int(wid))
                typ = form[-1].strip()
                dtype[label] = form_conv[typ]
                

        widths[1:-1] = [w+1 for w in widths[1:-1]]

        return headers, widths , dtype
        

    def read_state(self, state_filename):
        from .exomolstate import ExomolStates
        

        return ExomolStates(state_filename, self._pandas_state_fwf())

        
    @property
    def maximumTemperature(self):
        return self._max_temp

    @property
    def maximumPartitionTemperature(self):
        return self._max_temp_q

    @property
    def maximumWavenumber(self):
        return self._max_wavenumber
    
    @property
    def filePrefix(self):
        return f'{self._molecule_slug}__{self._linelist_name}'

    @property
    def availableBroadeners(self):
        return list(self._broadener_defs.keys())
    
    def create_broadeners(self, broadener):
        if broadener not in self.availableBroadeners:
            raise KeyError(f'Broadener with name {broadener} not available')
        else:
            return self._broadener_defs[broadener].generate_input() 


    def create_exocross_input(self, path='.'):
        from .exocrosswriter import ExocrossInput

        ex = ExocrossInput(self._molecule_slug, linelist=self._linelist_name, 
                           path=path, file_prefix=self.filePrefix)
        ex.set_molar_mass(self._mass)
        ex.set_range([0.01, self.maximumWavenumber])

        return ex
