

def run_pyexocross():
    import numpy as np
    import argparse
    from .pyexocross import PyExocross
    from .util import create_grid_res, convert_to_wavenumber
    parser = argparse.ArgumentParser()
    parser.add_argument("--linelist",type=str,dest="linelist",required=True)
    parser.add_argument("--path",type=str,dest="path",required=True)
    parser.add_argument("-T",type=float,dest="T",required=True)
    parser.add_argument("-P",type=float,required=True)
    parser.add_argument("--pressure-unit","-U",type=str,dest="U",default="bar",help="Units defined for pressure")
    parser.add_argument("-b","--broadeners",dest="broadeners",nargs="+",help='Broadeners to include')
    parser.add_argument("-r","--ratios",dest="ratios",nargs="+",type=float,help="corresponding ratios for each broadener. Default equally weighs them")
    parser.add_argument("-n","--nworkers",dest="nworkers",type=int,default=2,help="Number of worker threads to spin up for voigt calculation")
    parser.add_argument("-c","--chunk",type=int,default=100000,dest="chunk",help='How many transitions to read at a time')
    parser.add_argument("-o",type=str,dest="output",help="Output filename")
    parser.add_argument("-s",type=float,nargs="+",default=[0.1,10000], help='Spectral range')
    parser.add_argument("--thresh",type=float,default=1e-30, help='Threshold for intensities (default: %(default)s cm/molecule)')
    parser.add_argument("--wing",type=float,default=25.0, help='Voigt wing cutoff (default: %(default)s cm-1)')
    parser.add_argument("-u","--spectra-units",type=str,default="k",dest="spectra_unit",help="Spectral units, must be parsable by astropy (Default: cm-1)")
    parser.add_argument("-R","--resolution",dest="res",type=float,default=10000,help='Resolution (default: R=%(default)s)')
    parser.add_argument("--plot", action='store_true', default=False)

    args = parser.parse_args()

    linelist = args.linelist.lower()
    linelist_klass = None
    if linelist in ('exomol',):
        from .exomol import ExomolLinelist
        linelist_klass = ExomolLinelist
        print('Using ExoMol linelist')
    elif linelist in ('hitran',):
        from .hitran import HITRANLinelist
        linelist_klass = HITRANLinelist
        print('Using HITRAN linelist')
    else:
        raise ValueError(f'Unknown linelist {args.linelist}')
    
    ll = linelist_klass(args.path)

    print(f'Detected molecule is {ll.molecule}')


    temperature = args.T

    print(f'Temperature selected = {temperature} K')
    pressure_unit = args.U
    pressure_value = float(args.P)
    if pressure_unit != 'bar':
        from .util import conversion_factor
        factor = conversion_factor(pressure_unit,'bar')
        pressure_value*= factor
        print(f'Pressure selected {args.P} {pressure_unit} -> {pressure_value} bar')
    else:
        print(f'Pressure selected {pressure_value} bar')
    
    broadeners = args.broadeners
    ratios = args.ratios

    if broadeners is not None:
        if ratios is None:
            ratios = [1.0]*len(broadeners)
        if linelist in ('hitran',):
            for b,r in zip(broadeners,ratios):
                if b.lower() == 'self':
                    
                    ll.add_self_broadener(ratio=r)
                elif b.lower() =='air':
                    ll.add_air_broadener(ratio=r)
                else:
                    raise ValueError(f'HITRAN does not support broadener type {b}')
                print(f'Including {b} broadener at ratio={r}')
        if linelist in ('exomol',):
            for b,r in zip(broadeners,ratios):
                
                if b.lower() == 'default':
                    ll.add_default_broadener(ratio=r)
                elif b.lower() in ll.availableBroadeners:
                    ll.add_available_broadener(b, ratio=r)
                else:
                    raise ValueError(f'Exomol: Unknown {b} or file not found in linelist path')            
                print(f'Including {b} broadener at ratio={r}')

    max_jobs = args.nworkers*10

    R = args.res

    grid = convert_to_wavenumber(create_grid_res(R,min(args.s),max(args.s))[:,0],args.spectra_unit)
    grid = np.sort(grid)

    print(f'Running on wavenumber grid {grid.min()}--{grid.max()} cm-1 at R={R}')

    pyexo = PyExocross(ll)

    wn,xsec = pyexo.compute_xsec_parallel(grid,temperature,pressure_value, chunksize=args.chunk, threshold=args.thresh, wing_cutoff=args.wing,
                                          max_workers=args.nworkers,max_jobs=max_jobs)

    filename = args.output
    if filename is None:
        filename = f'{ll.molecule}_{temperature}K_{pressure_value}bar_R={R}.xsec'

    print(f'Writing output to {filename}')
    np.savetxt(filename, np.vstack((wn,xsec)).T)
    if args.plot:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(wn,xsec)
        plt.xlabel(r'Wavenumber cm$^{-1}$')
        plt.ylabel(r'Cross-section cm$^{2}$/molecule')
        plt.show()




