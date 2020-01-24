import numpy as np
import h5py

def get_filename(base_name,grid_x,grid_y,size,end_time,suffix,grid_z=None):
    if grid_z == None:
        return str(base_name) + "_ensemble=" + str(size) + "_" + str(grid_x) + "_" + str(grid_y) + "_" + str(end_time) + suffix + ".h5"
    else:
        return str(base_name) + "_ensemble=" + str(size) + "_" + str(grid_x) + "_" + str(grid_y) + "_" + str(grid_z) + "_" + str(end_time) + suffix + ".h5"

def get_path(directory, filename):
    return directory + filename

def py_out(pyfile,py_dataset,time):
    return pyfile[str(py_dataset)][str(py_dataset)+time][:]

#inner = (slice(2,-2),slice(2,-2))
inner = (slice(None,),slice(None,))

def ensemble_test_case(time, path, N, attribute, tag):

    file = h5py.File(path,'r')

    array = []
    for n in range(N):
        
        if type(time) == 'float':
            t_label = '_ensemble_mem=%i_%.3d_after_full_step' %(n,time)
        else:
#             t_label = '_ensemble_mem=%i_%s_after_full_step' %(n,time)
            t_label = '_%.3d_%s' %(time,tag)
        array.append(py_out(file,attribute,time=t_label)[inner])

    array = np.array(array)
    array = array.mean(axis=0)

    file.close()
    file.close()
    return array
