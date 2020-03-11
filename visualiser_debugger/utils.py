import numpy as np
import h5py

def cb_suffix(fs,ts):
    return "_cont_blend_fs=%i_ts=%i" %(fs,ts)

def get_tag_dict():
    td = {
            0 : 'before_flux',
            1 : 'before_advect',
            2 : 'after_advect',
            3 : 'after_ebnaexp',
            4 : 'after_ebnaimp',
            5 : 'after_half_step',
            6 : 'after_efna',
            7 : 'after_full_ebnaexp',
            8 : 'after_full_step',
    }
    return td
    
def spatially_averaged_rmse(arr,ref):
    arr = arr[2:-2,2:-2]
    ref = ref[2:-2,2:-2]
    
    arr -= arr.mean()
    ref -= ref.mean()

    return np.sqrt(((arr - ref)**2).mean())

def get_filename(base_fn,grid_x,grid_y,size,end_time,suffix):
    return base_fn + "_ensemble=" + str(size) + "_" + str(grid_x) + "_" + str(grid_y) + "_" + str(end_time) + suffix + ".h5"

def get_path(py_directory,filename):
    return py_directory + filename

def py_out(pyfile,py_dataset,time):
    return pyfile[str(py_dataset)][str(py_dataset)+time][:]

def test_case(time,path, N, attribute, label_type='TIME', tag='after_full_step', inner='None'):
    if inner == 'None':
        inner = (slice(None,),slice(None,))
    else:
        inner = (slice(2,-2),slice(2,-2))
        
    file = h5py.File(path,'r')
    
    if label_type == 'TIME':
        t_label = '_ensemble_mem=%i_%.3f_%s' %(n,time, tag)
    elif label_type == 'STEP':
        if N==1:
            t_label = '_%.3d_%s' %(time, tag)
        else:
            t_label = '_ensemble_mem=%i_%.3d_%s' %(n,time, tag)
            
    array = py_out(file,attribute,t_label)[inner]
    
    file.close()
    return np.array(array)

def ensemble_test_case(time, path, N, attribute, label_type='TIME', tag='after_full_step',inner='None'):
    if inner == 'None':
        inner = (slice(None,),slice(None,))
    else:
        inner = (slice(2,-2),slice(2,-2))
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
            
        array.append(py_out(file,attribute,time=t_label)[inner])

    array = np.array(array)
    array = array.mean(axis=0)

    file.close()
    file.close()
    return np.array(array)

def bin_func(obs,ens_mem_shape):
    obs = obs.reshape(ens_mem_shape[0],obs.shape[0]//ens_mem_shape[0],
                      ens_mem_shape[1],obs.shape[1]//ens_mem_shape[1])
    return obs.mean(axis=(1,3))

def rmse(diff):
    return np.sqrt((diff**2).mean())
