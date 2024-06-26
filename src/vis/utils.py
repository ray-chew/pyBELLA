import numpy as np
import h5py
import time

class test_case(object):
    def __init__(self,base_fn,py_dir,Nx,Ny,end_time,Nz=None):
        self.base_fn = base_fn
        self.py_dir = py_dir
        self.grid_x = Nx
        self.grid_y = Ny
        self.grid_z = Nz
        self.end_time = end_time
        if Nz is not None and Ny > 1:
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
        if self.grid_z is not None and Ny == 1:
            self.i2 = tuple([slice(2,-2)]*(self.ndim+1))
        else:
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
            'hcenter' : 'hcenter',
            'wplusx' : 'wplusx',
            'wplusy' : 'wplusy',
            'wplusz' : 'wplusz',
            'rhs' : 'rhs',
            'rhs_nodes' : 'rhs_nodes',
            'p2_full' : 'p2_full'
        }
        return dd

    def get_filename(self,N,suffix,format='h5'):
        if self.ndim == 2:
            fn = "%s_ensemble=%i_%i_%i_%.6f_%s.%s" %(self.base_fn,N,self.grid_x,self.grid_y,self.end_time,suffix,format)
        if self.ndim == 3 or self.grid_z is not None:
            fn = "%s_ensemble=%i_%i_%i_%i_%.6f_%s.%s" %(self.base_fn,N,self.grid_x,self.grid_y,self.grid_z,self.end_time,suffix,format)
        return fn


    def get_path(self,fn):
        path = self.py_dir + fn
        return path


    @staticmethod
    def py_out(pyfile,py_dataset,time):
        return pyfile[str(py_dataset)][str(py_dataset)+time][:], pyfile[str(py_dataset)][str(py_dataset)+time].attrs.get('t')


    def get_arr(self, path, time, N, attribute, label_type='TIME', tag='after_full_step', inner=False, avg=False, file=None):
        if inner == False:
            inner = self.i0
        else:
            inner = self.i2
            
        if file is None:
            file = h5py.File(path,'r')
        
        array = []
        
        if not hasattr(self,'t_arr'): self.t_arr = []
        t_arr = self.t_arr
        for n in range(N):
            if label_type == 'TIME':
                t_label = '_ensemble_mem=%i_%.3f_%s' %(n,time, tag)
            elif label_type == 'WINDOW_STEP':
                if N==1:
                    t_label = '_%.3d_%s' %(time, tag)
            elif label_type == 'STEP':
                if N==1:
                    t_label = '_ensemble_mem=0_%.3d_%s' %(time, tag)
                else:
                    t_label = '_ensemble_mem=%i_%.3d_%s' %(n,time, tag)
                
            arr, t = self.py_out(file,attribute,time=t_label)
            array.append(arr[inner])
            t_arr.append(t)

        array = np.array(array)
        if avg == True:
            array = array.mean(axis=0)

        if file is None:
            file.close()
        self.t_arr = t_arr
        return np.array(array)
    
    def spatially_averaged_rmse(self, arrs,refs,avg=False, grid_type='c'):
        diff = []
        refs = refs[:,np.newaxis,...]
        refs = np.repeat(refs, arrs.shape[1], axis=1)
#         print(arrs.shape, refs.shape)
        for arr, ref in zip(arrs,refs):
            #arr = (arr[self.i2])
            #ref = (ref[self.i2])
            if grid_type == 'n':
#                 print("arr before ie1 shape", arr.shape)
                arr = self.get_ie1(arr)
                ref = self.get_ie1(ref)
#                 print(arr.shape)
            if avg==True:
                for ens_mem in arr:
                     ens_mem -= ens_mem.mean()
#                 arr -= arr.mean()
                ref -= ref.mean()
                
            ref_ampl = ref.max() - ref.min()
            factor = ref_ampl
            factor = 1.0

            diff.append( np.sqrt(((arr - ref)**2).mean()) )
#             diff.append(np.sqrt( ((arr - ref)**2).mean() / (ref[0]**2).mean() ) )
#             diff.append(np.sqrt((((arr - ref) / ref)**2).mean()) )
            # Method A
            # err = [np.linalg.norm(mem - ref[0]) / np.linalg.norm(ref[0]) for mem in arr]
            #err = np.array(err).mean()
            # diff.append(err)
            # Method B
            # diff.append(np.linalg.norm(arr-ref) / np.linalg.norm(ref))
            # Method C
#             diff.append(np.linalg.norm((arr-ref).mean(axis=0)) / np.linalg.norm(ref[0]))
            
        return np.array(diff)

    def ensemble_spread(self, arrs,avg=False, grid_type='c'):
        diff = []
        
        refs = arrs.mean(axis=1)
        refs = refs[:,np.newaxis,...]
        refs = np.repeat(refs,arrs.shape[1],axis=1)
#         print(arrs.shape, refs.shape)

        for arr, ref in zip(arrs,refs):
            if grid_type == 'n':
                arr = self.get_ie1(arr)
                ref = self.get_ie1(ref)
            if avg==True:
                for ens_mem in arr:
                     ens_mem -= ens_mem.mean()
#                 arr -= arr.mean()
                ref -= ref.mean()
                
            diff.append(np.sqrt(((arr - ref)**2).mean()))
        return np.array(diff)
    
    
    def probe_rmse(self, arrs, refs, probe_loc, avg=False, inner=False):
        diff = []
       
        for arr, ref in zip(arrs,refs):
            if avg == True:
                arr = arr.mean(axis=0)
                ref = ref.mean(axis=0)

            if arr.ndim == 3:
                arr = arr[:,0,:]
                ref = ref[:,0,:]
               
            if inner == True:
                arr = arr[probe_loc[0],probe_loc[1]]
                ref = ref[probe_loc[0],probe_loc[1]]
            else:
                arr = arr[probe_loc[0],probe_loc[1]]
                ref = ref[probe_loc[0],probe_loc[1]]
                
            #diff.append(np.sqrt(((arr - ref)**2).mean()))
            diff.append(np.linalg.norm(arr-ref))
        return np.array(diff)
    
    # the first and last node rows / columns are repeated in periodic bcs, when we take mean we want to avoid this.
    @staticmethod
    def get_ie1(arr):
        arr = arr.squeeze()
        ndim = arr.ndim
        ie1 = tuple([slice(0,-1)]*ndim)
        arr = arr[ie1]
        return arr
    
    def get_mean(self, arrs, grid_type='c'):
        arr_mean = []
        
        for arr in arrs:
            if grid_type == 'n':
                arr = self.get_ie1(arr)
            arr_mean.append(arr.mean())
        
        return arr_mean
    

    @staticmethod
    def get_probe_loc(arrs, probe_loc):
        time_series = []
        for arr in arrs:
            time_series.append(arr[probe_loc[0],probe_loc[1]])
        return time_series
    
    
    def get_ensemble(self, times, N, attribute, suffix, cont_blend=False, ts=0, fs=0, label_type='TIME', tag='after_full_step', avg=False, diff=False, inner=True, load_ic=True, get_fn=True, fn=""):
        self.t_arr = []
        if cont_blend == True:
            suffix += cb_suffix(fs,ts)
            
        if get_fn:
            fn = self.get_filename(N,suffix)
        else:
            fn = fn
        path = self.get_path(fn)
        
        file = h5py.File(path,'r')
        
#         arr_lst = np.zeros((times.size+1),dtype=np.ndarray)
        arr_lst = []
        if load_ic:
            arr = self.get_arr(path, 0, N, attribute, tag='ic', label_type=label_type, avg=avg, inner=inner, file=file)
            arr_lst.append(arr)
        for tt, time in enumerate(times):
            arr = self.get_arr(path, time, N, attribute, tag=tag, label_type=label_type, avg=avg, inner=inner, file=file)
#             arr_lst[tt] = arr
            arr_lst.append(arr)
        
        
        if diff == True:    
            arr_lst = get_diff(arr_lst)
            
        file.close()    
        if load_ic:
            self.t_arr[0] = 0.0
        return np.array(arr_lst)
    
    
    def get_time_series(self, times, N, attribute, suffix, probe_loc, tag='after_full_step', cont_blend=False, ts=0, fs=0, label_type='TIME', diff=False, slc=[None], avg=False):
        if self.ndim == 3 and slc[0] == None:
            assert 0, "3D array has no 2D slice, define argument slc=(slice(...),slice(...),slice(...))"
            
        probe_row = probe_loc[0]
        probe_col = probe_loc[1]
        
        if cont_blend == True:
            suffix += '_' + self.cb_suffix(fs,ts)
            
        fn = self.get_filename(N,suffix)
        path = self.get_path(fn)
        
        probe = []
        for time in times:
            arr = self.get_arr(path, time, N, attribute, tag=tag, label_type=label_type, inner=True)
            
            if avg == True:
                arr -= arr.mean()
            
            if self.ndim == 3:
                arr = arr[slc].squeeze()
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
    # simple profiler for utils and plottting_tools
    def __init__(self, debug=True):
        self.tic = time.time()
        self.debug = debug
        
    def prtt(self, label=""):
        curr_time = time.time()
        if self.debug==True:
            print(label, curr_time - self.tic)
        self.tic = curr_time
