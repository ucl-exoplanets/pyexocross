import h5py
import datetime
import numpy as np
class HDF5Writer:


    def __init__(self, filename, molecule_name,T,P,append=False, flush_rate=20):
        self.filename = filename
        self.molecule_name = molecule_name
        self.fd = None
        self.P = []
        self.T = []
        self.temperatures = T
        self.pressures = P
        self.wngrid = None
        self.xsecs = []
        self.flushrate= flush_rate
        self.append = False
        self._write_count = 0
    def add_cross_section(self, T, P, wngrid, xsec):

        self.P.append(np.where(P == self.pressures)[0][0])
        self.T.append(np.where(T == self.temperatures)[0][0])
        self.xsecs.append(xsec)
        self._write_count += 1
        self.wngrid = wngrid
        if self._write_count > self.flushrate:
            self.flush_data()
            print('Flushing')
            self._write_count = 0

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, type, value, tb):
        self.close()

    def _openFile(self, fname):

        mode = 'w'

        fd = h5py.File(fname, mode=mode)
        fd.attrs['file_name'] = fname
        fd.attrs['file_time'] = datetime.datetime.now().isoformat()
        fd.attrs['creator'] = self.__class__.__name__
        fd.attrs['HDF5_Version'] = h5py.version.hdf5_version
        fd.attrs['h5py_version'] = h5py.version.version
        fd.attrs['program_name'] = 'pyExocross'
        fd.attrs['program_version'] = 'v3.0'
        return fd

    def open(self):
        self.fd = self._openFile(self.filename)
    

    def flush_data(self):



        if 'p' not in self.fd:
            self.fd.create_dataset('bin_edges', data=self.wngrid)
            self.fd.create_dataset('mol_name',data=self.molecule_name)
            self.fd.create_dataset('p',data=self.pressures)
            self.fd['p'].attrs['units'] = 'bar'
            self.fd.create_dataset('t',data=self.temperatures)
            
            self.fd.create_dataset('xsecarr',(len(self.pressures),len(self.temperatures),len(self.wngrid)))
        

        for P,T,xsecs in zip(self.P, self.T, self.xsecs):
            self.fd['xsecarr'][P,T] = xsecs

        self.fd.flush()
        self.P = []
        self.T = []
        self.xsecs = []




    def close(self):
        if self.fd:
            self.flush_data()
            self.fd.flush()
            self.fd.close()



