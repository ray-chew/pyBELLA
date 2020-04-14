import numpy as np
import h5py

class da_params(object):

    def __init__(self,N,da_type='rloc'):
        # number of ensemble members
        self.N = N
        # self.da_times = np.arange(0.0,3.75,0.25)[1:]
        self.da_times = np.arange(0.0,10.25,0.25)[1:]
        # self.da_times = np.arange(0.0,6.1,0.25)[1:]
        # self.da_times = []
        self.obs_attributes = ['rho','rhou','rhov','rhoY','p2_nodes']
        # self.obs_attributes = ['rho','rhou','rhov']
        # self.obs_attributes = ['rhoY','p2_nodeps']
        # self.obs_attributes = ['rhov']

        # self.obs_attributes = ['rhoY']
        # self.obs_attributes = ['rho', 'rhou', 'rhov','rhoY','p2_nodes']

        # which attributes to inflate in ensemble inflation?
        self.attributes = ['rho', 'rhou', 'rhov']

        # self.obs_path = './output_travelling_vortex/output_travelling_vortex_ensemble=1_32_32_6.0_truthgen.h5'
        self.obs_path = './output_rising_bubble/output_rising_bubble_ensemble=1_100_50_10.0_psinc_delth_ref.h5'

        # forward operator (projector from state space to observation space)
        self.forward_operator = np.eye(N)

        # localisation matrix
        # self.localisation_matrix = np.eye(N)
        weights = [0.05719096,0.25464401,0.33333333,0.25464401,0.05719096,0.25464401,0.52859548,0.66666667,0.52859548,0.25464401,0.33333333,0.66666667,1.,0.66666667,0.33333333,0.25464401,0.52859548,0.66666667,0.52859548,0.25464401,0.05719096,0.25464401,0.33333333,0.25464401,0.05719096]
        # weights = [0.01831564,0.082085,0.13533528,0.082085,0.01831564,0.082085,0.36787944,0.60653066,0.36787944,0.082085,0.13533528,0.60653066,1.,0.60653066,0.13533528,0.082085,0.36787944,0.60653066,0.36787944,0.082085,0.01831564,0.082085,0.13533528,0.082085,0.01831564]

        # weights3 = weights*len(self.obs_attributes)
        # self.localisation_matrix = np.diag(weights)
        # self.localisation_matrix += np.diag(np.ones_like(weights))
        # self.localisation_matrix = np.diag(np.ones((len(weights3))))
        self.localisation_matrix = weights
        self.localisation_matrix = np.ones_like(weights) * 1.0 #+ weights
        
        self.aprior_error_covar = 0.0001
        self.da_type = da_type

        # ensemble inflation factor for LETKF
        self.inflation_factor = 1.0

        # rejuvenation factor for ETPF
        self.rejuvenation_factor = 0.001

        self.loc = {
            'rho' : 0,
            'rhou' : 0,
            'rhov' : 0,
            'rhow' : 0,
            'rhoY' : 0,
            'rhoX' : 0,
            'p2_nodes' : 2,
        }

        self.loc_c = 0 # container list location of cell-based arrays
        self.loc_f = 1 # ... of face-based arrays
        self.loc_n = 2 # ... of node-based arrays

    def load_obs(self,obs_path,loc=0):
        if self.N > 1:
            obs_file = h5py.File(obs_path, 'r')
            obs_attributes = self.obs_attributes

            #### when were these observations taken?
            times = self.da_times
            # times = [1.0,2.0,3.0,4.0,5.0]

            #### axis 0 stores time series
            obs = np.empty(len(times), dtype=object)
            t_cnt = 0
            for t in times:
                #### how were these dataset called?
                label = '_ensemble_mem=0_%.3f_after_full_step' %t
                #### axis 1 stores the attributes
                obs[t_cnt] = {}
                for attribute in obs_attributes:
                    dict_attr = {
                        attribute: obs_file[str(attribute)][str(attribute) + str(label)][:]
                    }
                    obs[t_cnt].update(dict_attr)
                t_cnt += 1
            obs = np.array(obs)
            obs_file.close()
        return obs

    @staticmethod
    def sampler_gaussian(var):
        return lambda ic: np.random.normal(ic,var**0.5)

    @staticmethod
    def sampler_none():
        return lambda ic: ic

    @staticmethod
    def sampler_perturbator(max_shift):
        return lambda ic: np.roll(np.roll(ic, np.random.randint(-max_shift,max_shift), axis=1), np.random.randint(-max_shift,max_shift), axis=0)