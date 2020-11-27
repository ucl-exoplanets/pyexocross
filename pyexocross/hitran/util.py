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