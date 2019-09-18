import numpy as np
import h5py
import subprocess

class get_init_vals(object):
    def __init__(self):
        #####################################################
        # C-output folder names
        #
        self.test_cases = {2:'_acoustic_wave_high/',\
                            3:'_internal_long_wave/',\
                            1:'_travelling_vortex_3d_48_with_initial_projection/',\
                            4:'_output_rising_bubble/'
                            }

        self.base_suffix = self.test_cases[2]
        # 
        ####################################################

        #####################################################
        #
        # Time labels for C and Python
        #

        self.timestep = '000'
        self.t_label = 'after_full_step'
        # self.t_label_value = 9
        # self.labels = {0:self.t_label + '_' + 'ic',\
        #                 1:self.t_label + '_' + 'before_advect',\
        #                 2:self.t_label + '_' + 'after_advect',\
        #                 3:self.t_label + '_' + 'after_ebnaexp',\
        #                 4:self.t_label + '_' + 'after_ebnaimp',\
        #                 5:self.t_label + '_' + 'after_half_step',\
        #                 6:self.t_label + '_' + 'after_efna',\
        #                 7:self.t_label + '_' + 'after_full_advect',\
        #                 8:self.t_label + '_' + 'after_full_ebnaexp',\
        #                 9:self.t_label + '_' + 'after_full_step'
        #                 }
        self.label = self.timestep + '_' + self.t_label

        # self.label = self.labels[9]
        self.time = self.label
        # 
        ####################################################

    def update_test_case(self,test_case):
        self.base_suffix = self.test_cases[test_case]

    def update_time_step(self,timestep):
        self.timestep = timestep
        self.update_label(self.t_label)

    def update_label(self,label):
        self.label = self.timestep + '_' + label
        self.t_label = label

class converter(get_init_vals):
    def __init__(self):
        super().__init__()
        self.base_folder_name = "output"
        self.folder_comp = "low_Mach_gravity_comp/"
        self.folder_psinc = "low_Mach_gravity_psinc/"

        self.base_path = self.base_folder_name + self.base_suffix

        # specify the format extensions for input and output of converter
        self.hdf_format = '.'+ 'hdf'
        self.h5_format = '.' + 'h5'

        # Manually list the folder names for each output
        self.directories = ['S', 'T', 'Y', 'buoy', 'dT', 'dY', 'dp2_c', 'dp2_nodes', 'dpdim', 'drhoY', 'p', 'p2_c', 'p2_nodes', 'rho', 'rhoY', 'rhoe', 'u', 'v', 'vortz', 'w']
        self.directories += ['rhou', 'rhov', 'rhow', 'rhs']
        self.directories += ['buoy']

        # And manually list the file names.
        self.filenames = ['S', 'T', 'Y', 'buoy', 'dT', 'dY', 'dp2_c', 'dp2_n', 'dpdim', 'drhoY', 'p', 'p2_c', 'p2_n', 'rho', 'rhoY', 'rhoe', 'u', 'v', 'vortz', 'w']
        self.filenames += ['rhou', 'rhov', 'rhow', 'rhs']
        self.filenames += ['buoy']

        if self.time == self.t_label + '_' + 'after_ebnaimp':
            self.directories += ['wplusx', 'wplusy', 'hcenter', 'pnew', 'rhs_nodes', 'p2_initial']
            self.filenames += ['wplusx', 'wplusy', 'hcenter', 'p2_full', 'rhs_nodes', 'p2_initial']

        self.string_hdf = '_' + self.time + self.hdf_format
        self.string_h5 = '_' + self.time + self.h5_format
        self.full_path = self.base_path + self.folder_comp

        lst_hdf = []
        lst_h5 = []

        # build list of paths from parameters specified above.
        i = 0
        for directory in self.directories:
            lst_hdf.append(self.full_path + directory + '/' + self.filenames[i] + self.string_hdf)
            lst_h5.append(self.full_path + directory + '/' + self.filenames[i] + self.string_h5)
            i += 1 

        # print the list of paths
        # for path in lst_hdf:
        #     print(path)

        # now, convert the arrays specified in the list of paths
        # for item in lst_hdf:
        #     p = subprocess.call(["./h4toh5convert", item])