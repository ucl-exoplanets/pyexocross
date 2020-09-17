import subprocess
import tempfile
import os
import numpy as np
class ExocrossRunner:

    def __init__(self, exocross_path):
        self._exocross_path = exocross_path
        if not os.path.isfile(self._exocross_path):
            raise ValueError(f'Exocross path {self._exocross_path} does not exist')
    

    def run(self, exocross_input):

        with tempfile.TemporaryDirectory() as tmpdirname:
            output_filename = os.path.join(tmpdirname, 'output')
            input_filename = os.path.join(tmpdirname, 'input.par')

            stdin = exocross_input.generate_input(output_filename)
            print(stdin)
            subprocess.run([self._exocross_path], input=stdin,encoding='ascii')

            output = np.loadtxt(output_filename+'.xsec')

            return output


        
