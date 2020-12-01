
import matplotlib.pyplot as plt
from pyexocross import ExomolDef, ExocrossRunner, ExomolStates, ExomolTransitionReader
import numpy as np
import os
import time
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

states = exomol_def.read_state(exocross_input._state_file)

if os.path.isfile('test.csv'):
    os.unlink('test.csv')

total_transitions = 0

start = time.time()

final_v = []
final_transitions = []

for filename in exocross_input.valid_transitions():
    print(filename)
    for trans in ExomolTransitionReader.read_transitions(filename, chunksize=10000000):
        t = states.transition_states(trans, threshold=1e-28)
        total_transitions += len(t.index)
        t.to_csv('test.csv',index=False, header=False, mode='a', columns=['vif','Iif'])

total_time = time.time() - start

print(f'Total time for {total_transitions} transitions is {total_time} s. Average is {total_time/total_transitions} s')
