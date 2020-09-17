import matplotlib.pyplot as plt
from pyexocross import ExomolDef, ExocrossRunner
import numpy as np

path_to_exocross = '/Users/ahmed/Documents/repos/exocross/xcross.exe'

path_to_linelist = '/Users/ahmed/Documents/Linelists/H2O/'
path_to_def_file = '/Users/ahmed/Documents/Linelists/H2O/1H2-16O__BT2.def'


# Create our exocross runner

exo_run = ExocrossRunner(path_to_exocross)

# Read exomol-def file

exomol_def = ExomolDef(path_to_def_file)

print(f'Available Broadeners: {exomol_def.availableBroadeners}')

# Create an exocross input

exocross_input = exomol_def.create_exocross_input(path_to_linelist)

# Generate broadeners
h2_broad = exomol_def.create_broadeners('H2')
he_broad = exomol_def.create_broadeners('He')

# Add broadeners to input
exocross_input.add_broadener(h2_broad)
exocross_input.add_broadener(he_broad)

# Set their ratios
h2_broad.ratio = 0.9
he_broad.ratio = 0.1

# Set wavelength range and point
exocross_input.set_range([0.1, 2], units='um')
exocross_input.Npoints = 101

# Lets set some temperature and pressures to run
temperature_points = np.linspace(500, 2000, 10)
pressure_points = np.logspace(0, 6, 10)

final_grid = np.zeros(shape=(10,10,101))
wn_grid = None
for ip, p in enumerate(pressure_points):
    for it, t in enumerate(temperature_points):
        exocross_input.temperature = t
        exocross_input.pressure = p

        res = exo_run.run(exocross_input)
        if wn_grid is None:
            wngrid = res[:, 0]

        final_grid[ip, it] = res[:, 1]

plt.figure()
plt.plot(10000/wngrid, final_grid[-1, -1])
plt.show()
