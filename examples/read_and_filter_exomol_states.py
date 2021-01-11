import matplotlib.pyplot as plt
from pyexocross.exomol import ExomolDef
from pyexocross.exomol import ExomolLinelist

# Load an exomol def file and get all info

exomol_def = ExomolDef('/Users/ahmed/Documents/Linelists/H2O/1H2-16O__BT2.def')

# Use the exomol def to read in the state file. Can be either the compressed or uncompressed

# states = exomol_def.read_state('/Users/ahmed/Documents/Linelists/H2O/1H2-16O__BT2.states.bz2')

states = exomol_def.read_state('/Users/ahmed/Documents/Linelists/H2O/1H2-16O__BT2.states')

# Now we can get the pandas dataframe

df = states.pandasFrame

# We can print only states we care about e.g fitler for J > 10 and E > 3000 cm-1

print(df[ (df['J'] > 10) & (df['E'] > 3000) ]) 


# Now lets try using the entire linelist
# The linelist will automatically handle finding the def, states, pf, broads and transitions

linelist = ExomolLinelist(path_to_linelist='/Users/ahmed/Documents/Linelists/H2O/')

# Lets get our temperature
T = 1000

# Our wavenumber range to read from
min_wn, max_wn = 1000, 5000

# We will read 10000 transitions at a time
for tr in linelist.get_transitions(min_wn, max_wn, chunksize=10000):

    # Now we can filter the transitions for specific quanta. We get " and ' for initial and final respectively
    # Lets filter for all transitions J' > 10 and v1' == 1 and v" == 0
    filtered = tr[ (tr["J'"] > 10) & (tr["v1'"] == 1) & (tr['v1"'] == 0) ]

    # Compute their intensity
    I = linelist.compute_intensity(filtered, T)

    # Add the intensity to the transition

    filtered.insert(1,"I_if", I)
    
    # Filter low intensity
    filtered = filtered[ filtered['I_if'] > 1e-28]


    # Print or save the result to whatever: Excel, CSV etc!
    print(filtered)





