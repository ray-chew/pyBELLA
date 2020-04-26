import numpy as np
import h5py
from time import time

class test_case(object):
    def __init__(self,base_fn,py_dir,Nx,Ny,end_time,Nz=None):
        self.base_fn = base_fn
        self.py_dir = py_dir
        self.grid_x = Nx
        self.grid_y = Ny
        self.end_time = end_time
        if Nz != None:
            self.grid_z = Nz
            self.ndim = 3
        else:
            self.ndim = 2
        
        self.cb_suffix = self.cb_suffix
        self.get_tag_dict = self.get_tag_dict
        self.py_out = self.py_out
        
        self.get_filename = self.get_filename
        self.get_path = self.get_path
        self.get_arr = self.get_arr
        
        self.i0 = tuple([slice(None,)]*self.ndim)
        self.i2 = tuple([slice(2,-2)]*self.ndim)
        
    def cb_suffix(self,fs,ts,suffix=""):
        if suffix != "":
            return "%s_cont_blend_fs=%i_ts=%i" %(suffix,fs,ts)
        else:
            return "cont_blend_fs=%i_ts=%i" %(fs,ts)


    @staticmethod
    def get_tag_dict():
        td = {
                0 : 'before_flux',
                1 : 'before_advect',
                2 : 'after_advect',
                3 : 'after_ebnaexp',
                4 : 'after_ebnaimp',
                5 : 'after_half_step',
                6 : 'after_efna',
                7 : 'after_full_advect',
                8 : 'after_full_ebnaexp',
                9 : 'after_full_step',
        }
        return td

    @staticmethod
    def get_debug_attrs():
        dd = {
            'p2_initial' : 'p2_initial',
            'hecenter' : 'hcenter',
            'wplusx' : 'wplusx',
            'wplusy' : 'wplusy',
            'wplusz' : 'wplusz',
            'rhs' : 'rhs',
            'rhs_nodes' : 'rhs_nodes',
            'p2_full' : 'p2_full'
        }
        return dd

    def get_filename(self,N,suffix):
        if self.ndim == 2:
            fn = "%s_ensemble=%i_%i_%i_%.1f_%s.h5" %(self.base_fn,N,self.grid_x,self.grid_y,self.end_time,suffix)
        if self.ndim == 3:
            fn = "%s_ensemble=%i_%i_%i_%i_%.1f_%s.h5" %(self.base_fn,N,self.grid_x,self.grid_y,self.grid_z,self.end_time,suffix)
        return fn


    def get_path(self,fn):
        path = self.py_dir + fn
        return path


    @staticmethod
    def py_out(pyfile,py_dataset,time):
        return pyfile[str(py_dataset)][str(py_dataset)+time][:]


    def get_arr(self, path, time, N, attribute, label_type='TIME', tag='after_full_step', inner=False, avg=True):
        if inner == False:
            inner = self.i0
        else:
            inner = self.i2
            
        file = h5py.File(path,'r')
        
        array = []
        for n in range(N):
            if label_type == 'TIME':
                t_label = '_ensemble_mem=%i_%.3f_%s' %(n,time, tag)
            elif label_type == 'STEP':
                if N==1:
                    t_label = '_%.3d_%s' %(time, tag)
                else:
                    t_label = '_ensemble_mem=%i_%.3d_%s' %(n,time, tag)
                
            array.append(self.py_out(file,attribute,time=t_label)[inner])

        array = np.array(array)
        if avg == True:
            array = array.mean(axis=0)

        file.close()
        return np.array(array)
    
    @staticmethod
    def spatially_averaged_rmse(arrs,refs,avg=False):
        diff = []
        for arr, ref in zip(arrs,refs):
            arr = arr[self.i2]
            ref = ref[self.i2]
            
            if avg==True:
                arr -= arr.mean()
                ref -= ref.mean()

            diff.append(np.sqrt(((arr - ref)**2).mean()))
        return np.array(diff)
    
    @staticmethod
    def get_probe_loc(arrs, probe_loc):
        time_series = []
        for arr in arrs:
            time_series.append(arr[probe_loc[0],probe_loc[1]])
        return time_series
    
    
    def get_ensemble(self, times, N, attribute, suffix, cont_blend=False, ts=0, fs=0, label_type='TIME'):
        if cont_blend == True:
            suffix += cb_suffix(fs,ts)
            
        fn = self.get_filename(N,suffix)
        path = self.get_path(fn)
        
        arr_lst = []
        for time in times:
            arr = self.get_arr(path, time, N, attribute, label_type=label_type)
            arr_lst.append(arr)
            
        return np.array(arr_lst)
    
    
    def get_time_series(self, times, N, attribute, suffix, probe_loc, cont_blend=False, ts=0, fs=0, label_type='TIME', diff=False):
        probe_row = probe_loc[0]
        probe_col = probe_loc[1]
        
        if cont_blend == True:
            suffix += cb_suffix(fs,ts)
            
        fn = self.get_filename(N,suffix)
        path = self.get_path(fn)
        
        probe = []
        for time in times:
            arr = self.get_arr(path, time, N, attribute, label_type=label_type)
            probe.append(arr[probe_row,probe_col])
            
        probe = np.array(probe)
        if diff==True:
            return get_diff(probe)
        elif diff==False:
            return probe



def bin_func(obs,ens_mem_shape):
    obs = obs.reshape(ens_mem_shape[0],obs.shape[0]//ens_mem_shape[0],
                      ens_mem_shape[1],obs.shape[1]//ens_mem_shape[1])
    return obs.mean(axis=(1,3))

def rmse(diff):
    return np.sqrt((diff**2).mean())

def get_diff(probe):
    probe = np.array(probe)
    return probe[1:] - probe[:-1]

def spatially_averaged_rmse(arr,ref):
    #arr = arr[2:-2,2:-2]
    #ref = ref[2:-2,2:-2]
    
    #arr -= arr.mean()
    #ref -= ref.mean()

    return np.sqrt(((arr - ref)**2).mean())

class prt_time(object):
    def __init__(self, debug=True):
        self.tic = time()
        self.debug = debug
        
    def prtt(self, label=""):
        curr_time = time()
        if self.debug==True:
            print(label, curr_time - self.tic)
        self.tic = curr_time
