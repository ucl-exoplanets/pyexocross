def read_hitran_line(line):
    import numpy as np
    lengths = [2, 1, 12, 10, 10, 5, 5, 10, 4, 8, 15, 15]
    position = np.array([0] + lengths).cumsum()[:-1]

    formats = [int,int,float,float,float,float,float,float,float,float,str,str]
    res = [f(line[ls:ls+s]) for ls,s,f in zip(position, lengths, formats)] + \
        [float(line[-14:-7]), float(line[-7:])]
    return res


def read_hitran(filename):
    with open(filename,'r') as f:
        for line in f:
            yield read_hitran_line(line)


def read_hitran_pandas(hitran_file, chunksize=None):
    import pandas
    lengths = [2, 1, 12, 10, 10, 5, 5, 10, 4, 8, 15, 15, 15, 15, 6, 12, 1, 7, 7]
    names = ['MolID', 'IsoID', 'v_if','S_if','A_if','gamma_air','gamma_self','E"', 
             'n_air','delta_air',"g_quanta'",'g_quanta"',"l_quanta'",'l_quanta"','error','Refs',
             'Flag',"g'",'g"']
    df = pandas.read_fwf(hitran_file, widths=lengths, names=names, chunksize=2)
    chunk = next(df)
    molId = chunk['MolID'][0]
    quanta_names, quanta_widths = ['g_quanta',], [15,]
    res = moleculeClass(molId)
    quanta_names, quanta_widths, quanta_start,skip = res
    quanta_widths[0] += quanta_start

    lquanta_up, lwidths_up, lquanta_low, lwidths_low = moleculeGroup(molId)
    
    lengths = [2, 1, 12, 10, 10, 5, 5, 10, 4, 8, *quanta_widths, *quanta_widths, 
               *lwidths_up, *lwidths_low, 6, 12, 1, 7, 7]
    names = ['MolID', 'IsoID', 'v_if','S_if','A_if','gamma_air','gamma_self','E"', 
             'n_air','delta_air',*[f"{x}'" for x in quanta_names],
             *[f'{x}"' for x in quanta_names],*lquanta_up,*lquanta_low,'error','Refs',
             'Flag',"g'",'g"']
    
    df = pandas.read_fwf(hitran_file, widths=lengths, names=names, chunksize=chunksize)
    for chunk in df:
        yield chunk


def moleculeGroup(molid):
    from .hapi import moleculeName
    molname = moleculeName(molid)
    if molname in ('H2O','O3','SO2','NO2','HNO3','H2CO','HOCl','H2O2','COF2',
                   'H2S','HO2','HCOOH','ClONO2','HOBr','C2H4'):
        quanta = ['J','Ka','Kc','F','Sym']
        widths = [3, 3, 3, 5, 1]
        return [f"{x}'" for q in quanta], widths, \
            [f'{x}"' for x in quanta], widths

    if molname in ('CO2', 'N2O', 'CO', 'HF','HCl', 'HBr',
                   'HI','OCS','N2','HCN','C2H2','NO+'):
        upper_quanta = ["F'",]
        upper_widths = [15,]
        lower_quanta = ['Br','J"','Sym"','F"',]
        lower_widths = [6,3,1,5,]
        return upper_quanta, upper_widths, \
                lower_quanta,lower_widths
    if molname in ('SF6', 'CH4',):
        quanta = ['J','C','a','F',]
        widths = [5,2,3,5,]
        return [f"{x}'" for q in quanta], widths, \
            [f'{x}"' for x in quanta], widths
    if molname in ('CH3D','CH3Cl','C2H6', 'NH3','PH3',):
        quanta = ['J','Ka','l','C','Sym','F']
        widths = [3,3,2,2,1,4,]
        return [f"{x}'" for q in quanta], widths, \
            [f'{x}"' for x in quanta], widths
    if molename in ('O2',):
        upper_quanta = ["F'",]
        upper_widths = [15,]
        lower_quanta = ['Br1','N"','Br2','Sym"','J"','F"',]
        lower_widths = [2,3,1,3,1,5,]
        return upper_quanta, upper_widths, \
                lower_quanta,lower_widths  
    if molname in ('NO','OH','ClO'):
        upper_quanta = ["F'",]
        upper_widths = [15,]
        lower_quanta = ['Br', 'J"','Sym"','F"',]
        lower_widths = [4, 5, 1, 5,]
        return upper_quanta, upper_widths, \
                lower_quanta,lower_widths  
    
    return "l_quanta'",15,'l_quanta"',15
def moleculeClass(molid):
    from .hapi import moleculeName
    molname = moleculeName(molid)
    if molname in ('CO', 'HF', 'HCl', 'HBr', 'HI', 'N2','NO+'):
        return ['v1',],[2,],13, None
    if molname in ('O2',):
        return ['X', 'v1', ],[1, 2,],12, None
    if molname in ('NO','OH','ClO',):
        return ['X', 'i', 'v1', ],[1, 3, 4],7, None
    if molname in ('N2O', 'OCS', 'HCN',):
        return ['v1', 'v2','l2','v3',], [2,2,2,2], 7, None
    if molname in ('CO2',):
        return ['v1', 'v2','l2','v3','r', ], [2,2,2,2,1], 6, None

    if molname in ('H2O', 'O3', 'SO2', 'NO2','HOCl','H2S',
                   'HO2','HOBr',):
        return ['v1', 'v2', 'v3',],[2,2,2,],9,None
    if molname in ('C2H2',):
        return ['v1', 'v2', 'v3', 'v4', 'v5','l','+-','r', 'S',],[2,2,2,2,2,2,1,1,1],0,None
    
    if molname in ('NH3', 'PH3',):
        return ['v1','v2', 'v3', 'v4', 'S'],[2,2,2,2,2],5,None
    
    if molname in ('H2CO', 'H2O2','COF2',):
        return ['v1', 'v2', 'v3', 'v4', 'v5','v6'],[2,2,2,2,2,2], 3, None
    
    if molname in ('CH4', 'CH3D','CH3Cl','C2H6','HNO3',
                   'SF6','HCOOH','ClONO2','C2H4',):
        return ['v1', 'v2', 'v3', 'v4','n','C'],[2,2,2,2,2,2],3,None

    return ['g_quanta',],[15,],0,None


def splitwidths(string, widths, start=0, skip=None):
    import numpy as np
    position = np.cumsum([start] + widths[:-1])
    if skip is None:
        skip=[False]*len(widths)
    return [string[p:p+w] for p,w,s in zip(position,widths,skip) if not s]






