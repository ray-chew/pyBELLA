import numpy as np
import utils
import plotting_tools as pt
import h5py
import os

from scipy.interpolate import griddata

class load_gen(object):
    def __init__(self, Nx, Ny, et, N, Nz = None):
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
            
        self.et = et
        self.N = N
       
    def get_tc(self,base_fn, pydir = None):
        self.base_fn = base_fn
        if pydir is not None:
            pydir = pydir
        else:
            pydir = '../%s/' %base_fn

        tc = utils.test_case(base_fn, pydir, self.Nx, self.Ny, self.et, Nz=self.Nz)
        tags = tc.get_tag_dict()
        
        self.tc, self.tags = tc, tags
        return tc, tags
    
    def get_path(self, sfx):
        fn = self.tc.get_filename(self.N,sfx)
        path = self.tc.get_path(fn)
        return path
    
    def load_data(self, sfx, time_series, attributes, t_lbl = 'TIME', tag = None, inner = True, slc = 0):
        """
        Load data for truth generation. Currently works for 2D horizontal and vertical slices. 3D horizontal slices are projected down to 2D for conversion / truth generation.
        
        Returns
        -------
        ts : nd.array
            An array of size [(time series length) x (attributes length)]
        """
        
        time_series = np.array(time_series)
        self.time_series = time_series
        attributes = np.array(attributes)
        self.attributes = attributes
        
        if inner == False and slc == 0:
            slc = 2
        
        if self.Nz is not None:
            arr_slc = (slice(None,),slc,slice(None,))
        else:
            arr_slc = (slice(None,),slice(None,))
            
        path = self.get_path(sfx)
            
        ts = np.zeros((time_series.shape[0], attributes.shape[0]), dtype='object')
        
        if tag == None:
            tag = self.tags[9]
        
        for tt, time in enumerate(time_series):
            for aa, attribute in enumerate(attributes):
                ts[tt][aa] = self.tc.get_arr(path, time, self.N, attribute, label_type=t_lbl, tag=tag, inner=inner, avg=False)[0][arr_slc]
        
        self.ts = ts
        return ts
            
    def convert(self, method, nNx, nNy):
        """
        Convert does the generation of the truth from given data. Accepts method function / class that does the conversion.
        
        Note
        ----
        Works only in 2D for now.
        
        """
        None
        
    def debug(self):
        None
        
    def save(self, tsn, sfx):
        if self.Nz is None:
            file = h5py.File('../%s/%s_ensemble=%i_%i_%i_%.1f_%s.h5' %(self.base_fn, self.base_fn, self.N, self.Nx, self.Ny, self.et, str(sfx)))
        else:
            
            fn = '../%s/%s_ensemble=%i_%i_%i_%i_%.1f_%s.h5' %(self.base_fn, self.base_fn, self.N, self.Nx, self.Ny, self.Nz, self.et, str(sfx))
            if os.path.exists(fn):
                os.rename(fn, fn + '_old')
            
            file = h5py.File(fn, 'w')
            
        for tt, ts_data in enumerate(tsn):
            time = self.time_series[tt]
            name = "ensemble_mem=0_%.3f_after_full_step" %time
            
            for aa, attribute in enumerate(self.attributes):
                ts_attrdata = ts_data[aa]
                
                # Since the load method projects 3D horizontal slices -> 2D, if original data is 3D, recreate the vertical axis.
                if self.Nz is not None:
                    ts_attrdata = ts_attrdata[:, np.newaxis, :]
                    
                file.create_dataset(str(attribute) + '/' + str(attribute) + '_' + str(name), data=ts_attrdata, chunks=True, compression='gzip', compression_opts=4, dtype=np.float32)
        file.close()
        
       
class interpolate(object):
    def __init__(self, nNx, nNy, Lx, Ly, a_type):
        self.nNx = nNx
        self.nNy = nNy
        
        self.Lx = Lx
        self.Ly = Ly
        
        self.a_type = a_type
        
        # cells and nodes meshgrids for the new grid onto which data will be interpolated.
        self.gridc_x, self.gridc_y = np.meshgrid(np.linspace(0,Lx,nNx), np.linspace(0,Ly,nNy))
        self.gridn_x, self.gridn_y = np.meshgrid(np.linspace(0,Lx,nNx+1),np.linspace(0,Ly,nNy+1))
        
    def convert(self, ts, data_type='array'):
        """
        conversion / interpolation takes place here. 
        """
        
        if data_type == 'array':
            tlen = ts.shape[0] # length of time-series
            alen = ts.shape[1] # length of attributes
        
        tsn = np.zeros_like(ts)
        
        for aa in range(alen):
            if self.a_type[aa] == 'cell':
                grid_x, grid_y = self.gridc_x, self.gridc_y
                nx, ny = self.nNx, self.nNy
            elif self.a_type[aa] == 'node':
                grid_x, grid_y = self.gridn_x, self.gridn_y
                nx, ny = self.nNx + 1, self.nNy + 1
                
            for tt in range(tlen):
                arr = ts[tt][aa]
                
                gNx, gNz = arr.shape[0], arr.shape[1]
                x, z = np.linspace(0,self.Lx,gNx), np.linspace(0,self.Ly,gNz)
                X,Z = np.meshgrid(x,z)
            
                points = np.zeros((np.zeros((gNx,gNz)).flatten().shape[0],2))
                points[:,0] = X[...].flatten()
                points[:,1] = Z[...].flatten()
            
                values = (arr).flatten()

                arr_interp = griddata(points, values, (grid_x, grid_y), method='cubic')
                arr_interp = arr_interp.reshape(nx,ny)
                
                tsn[tt][aa] = arr_interp
                
        return tsn
