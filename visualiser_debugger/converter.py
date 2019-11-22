import subprocess
import numpy as np
import h5py

class convert(object):
    # specify the format extensions for input and output of converter
    hdf_format = '.'+ 'hdf'
    h5_format = '.' + 'h5'
    
    def __init__(self,base_path, folder_path, time, t_label):
        self.base_path = base_path
        self.folder_path = folder_path
        
        self.string_hdf = '_' + time + self.hdf_format
        self.string_h5 = '_' + time + self.h5_format
        self.full_path = base_path + folder_path
        
        self.t_label = t_label
        self.time = time
        
    def convert_files(self,directories,filenames):
        if self.time == self.t_label + '_' + 'after_ebnaimp':
            directories += ['wplusx', 'wplusy', 'hcenter', 'pnew', 'rhs_nodes', 'p2_initial']
            filenames += ['wplusx', 'wplusy', 'hcenter', 'p2_full', 'rhs_nodes', 'p2_initial']

        lst_hdf = []
        lst_h5 = []

        # build list of paths from parameters specified above.
        i = 0
        for directory in directories:
            lst_hdf.append(self.full_path + directory + '/' + filenames[i] + self.string_hdf)
            lst_h5.append(self.full_path + directory + '/' + filenames[i] + self.string_h5)
            i += 1 

        # print the list of paths
        for path in lst_hdf:
            print(path)

        # now, convert the arrays specified in the list of paths
        for item in lst_hdf:
            p = subprocess.call(["./h4toh5convert", item])
        
    def get_converted_files(self,case_folders,case_names):
        test_cases_folders = case_folders
        test_cases_names = case_names
        
        if self.time == self.t_label + '_' + 'after_ebnaimp':
            test_cases_folders += ['wplusx', 'wplusy', 'hcenter', 'pnew', 'rhs_nodes', 'p2_initial']
            test_cases_names += ['wplusx', 'wplusy', 'hcenter', 'p2_full', 'rhs_nodes', 'p2_initial']
            
        test_cases_folders = np.char.array(test_cases_folders)
        test_cases_names = np.char.array(test_cases_names)
            
        # build paths from folder names and filenames
        test_cases_paths = self.full_path + test_cases_folders + '/' + test_cases_names + self.string_h5

        # define empty class as a holder for all the C-hdf5 output.
        class c_output(object):
            def __init__(self):
                None
                
            def c_out(self):
                for key,value in vars(self).items():
                    setattr(self,key,value['Data-Set-2'])
                    
        # get an instance of an empty class to populate it with the C-HDF output as attributes
        c = c_output()

        # populate the class with the C-HDF output as attributes
        i = 0
        for path in test_cases_paths:
            setattr(c,test_cases_folders[i],h5py.File(path, 'r'))
            i += 1
            
        # run the method to extract array from HDF objects
        c.c_out()
        return c
