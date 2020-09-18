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
        val = bool(self._lines[self._count])
        self.skip(skip)
        return val

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

    def add_code(self, quanta_code):
        self._avail_codes.append(quanta_code)

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
        lr.skip(num_cases)
        no_quanta = lr.read_int()
        lr.skip(no_quanta * 3)
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
                    new_broad.add_code(code_label)
                    no_quanta = lr.read_int()
                    lr.skip(no_quanta)

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
