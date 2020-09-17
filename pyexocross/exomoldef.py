from taurex.log import Logger

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
        
        count = 0
        if self.exocross_lines[count] != 'EXOMOL.def':
            raise IOError('Incorrect EXOMOL def header')
        count += 2
        self._molecule_slug = self.exocross_lines[count]
        count += 1
        self._linelist_name = self.exocross_lines[count]
        count += 1
        self._version_number = self.exocross_lines[count]
        count += 1
        self._inchikey = self.exocross_lines[count]

        self.info(f'Molecule is {self._molecule_slug}')
        self.info(f'Linelist: {self._linelist_name} '
                  f'Version: {self._version_number}')
        count += 1
        while(True):
            line = self.exocross_lines[count]
            try:
                split = line.split()
                self._mass = float(split[0])
                float(split[1])
                break
            except (IndexError, ValueError, ):
                count += 1
                continue

        count += 2  # Skip symmetry
        num_irr = int(self.exocross_lines[count])
        count += num_irr * 3 + 1
        self._max_temp = float(self.exocross_lines[count])
        count += 1
        self._num_broadeners = int(self.exocross_lines[count])
        count += 6
        self._num_states = int(self.exocross_lines[count])
        count += 1
        num_cases = int(self.exocross_lines[count])
        count += num_cases + 1
        no_quanta = int(self.exocross_lines[count])
        count += no_quanta * 3 + 1
        self._total_transitions = int(self.exocross_lines[count])
        count += 1 
        self._num_trans_files = int(self.exocross_lines[count])
        count += 1
        self._max_wavenumber = float(self.exocross_lines[count])
        count += 1
        self._highest_complete = float(self.exocross_lines[count])
        count += 1
        self._max_temp_q = float(self.exocross_lines[count])
        count += 1
        self._t_step = float(self.exocross_lines[count])
        count += 2
        self._default_gamma = float(self.exocross_lines[count])
        count += 1
        self._default_n = float(self.exocross_lines[count])
        self._broadener_defs = {}

        if self._num_broadeners > 0:
            count += 1
            for b in range(self._num_broadeners):

                broadener_label = self.exocross_lines[count]
                self.info(f'Reading broadener {broadener_label}')
                count += 1
                broadener_filename = self.exocross_lines[count]
                count += 1
                jmax = int(self.exocross_lines[count])
                count += 1
                default_gamma = float(self.exocross_lines[count])
                count += 1
                default_n = float(self.exocross_lines[count])
                count += 1

                new_broad = BroadenerData(broadener_label, broadener_filename,
                                          jmax, default_gamma, default_n)
                self._broadener_defs[broadener_label] = new_broad

                n_broad_quanta_set = int(self.exocross_lines[count])
                count += 1
                for x in range(n_broad_quanta_set):
                    code_label = self.exocross_lines[count]
                    new_broad.add_code(code_label)
                    count += 2
                    no_quanta = int(self.exocross_lines[count])
                    count += no_quanta + 1

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
