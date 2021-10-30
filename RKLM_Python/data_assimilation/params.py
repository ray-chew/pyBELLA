import numpy as np
import h5py
from scipy import stats # generate Gaussian localisation matrix

from scipy import signal
from inputs.boundary import set_explicit_boundary_data, set_ghostnodes_p2 # for converter

class da_params(object):

    def __init__(self,N,da_type='rloc'):
        # number of ensemble members
        self.N = N

        self._da_times = np.arange(0.0,3.25,0.25)[1:]
        # self._da_times = np.arange(5.0,10.5,0.5)/10.0
        # self._da_times = [0.1]
        self._da_times = np.around(self.da_times,3)
        
        self.obs_attributes = ['rho','rhou', 'rhov', 'rhoY', 'p2_nodes']
        # self.obs_attributes = ['rhou','rhov']

        # which attributes to inflate in ensemble inflation?
        self.attributes = ['rho', 'rhou', 'rhov']

        # self.obs_path = './output_travelling_vortex/output_travelling_vortex_ensemble=1_32_32_6.0_truthgen.h5'
        # self.obs_path = './output_rising_bubble/output_rising_bubble_ensemble=1_100_50_10.0_psinc_ref.h5'
        # self.obs_path = './output_rising_bubble/output_rising_bubble_ensemble=1_100_50_10.0_truthgen_freezelt5.h5'
        # self.obs_path = './output_rising_bubble/output_rising_bubble_ensemble=1_100_50_10.0_comp_ref.h5'
        # self.obs_path = './output_rising_bubble/output_rising_bubble_ensemble=1_160_80_1.0_truth_CFLfixed_ib-0.h5'

        # self.obs_path = './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pps_tra_truth.h5'
        # self.obs_path = './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_neg_comp_1.0_pp_tra_truth_ip.h5'
        # self.obs_path = './output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_comp_1.0_pp_tra_truth_ip.h5'
        self.obs_path = './output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5'

        # forward operator (projector from state space to observation space)
        self.forward_operator = np.eye(N)
        # self.converter = self.converter

        ############################################
        # Parameters for sparse observations
        ############################################
        self.sparse_obs = True
        self.sparse_obs_by_attr = False

        if self.sparse_obs_by_attr:
            assert(0, "Not yet implemented.")

        self.obs_frac = 0.90 # fraction of the observations to pick.
        self.gen_obs_sparse()

        ############################################
        # Parameters for measurement noise
        ############################################
        self.add_obs_noise = True
        self.noise_type = 'VarCov'

        self._noise_percentage = 0.05
        self.obs_noise = {}
        self.std_dev = []
        self.gen_obs_noise()

        ############################################
        # Parameters for LETKF subdomain size
        ############################################
        self.obs_X = 11
        self.obs_Y = 11

        # constants, linear, gaussian
        self.localisation_matrix = self.get_loc_mat('gaussian')

        self.da_type = da_type

        # ensemble inflation factor for LETKF
        self.inflation_factor = 1.00

        # rejuvenation factor for ETPF
        self.rejuvenation_factor = 0.001

        self.loc_c = 0 # container list location of cell-based arrays
        self.loc_f = 1 # ... of face-based arrays
        self.loc_n = 2 # ... of node-based arrays

        # in which data container are the attributes involved in the DA procedure?
        self.loc = {
            'rho' : 0,
            'rhou' : 0,
            'rhov' : 0,
            'rhow' : 0,
            'rhoY' : 0,
            'rhoX' : 0,
            'p2_nodes' : 2,
        }

    def gen_obs_sparse(self):
        da_len = len(self.da_times)
        if self.sparse_obs_by_attr == True:
            da_depth = len(self.obs_attributes)
        else:
            da_depth = 1
        if self.sparse_obs == True and da_len > 0:
            np.random.seed(777)
            self.sparse_obs_seeds = np.random.randint(10000, size=(da_len, da_depth))
            if len(self.sparse_obs_seeds) > 1:
                self.sparse_obs_seeds = self.sparse_obs_seeds.squeeze()


    def gen_obs_noise(self):
        for cnt, key in enumerate(self.obs_attributes):
            if self.noise_type == 'FixCov':
                assert self.std_dev is not None, "std_dev keyword argument must be a list equal in size to dap.obs_attributes"
                assert len(self.std_dev) == len(self.obs_attributes), "std_dev keyword argument must be a list equal in size to dap.obs_attributes"
                self.obs_noise[key] = float(self.std_dev[cnt])
                print(self.std_dev[cnt])
                cnt += 1
            else:    
                self.obs_noise[key] = self._noise_percentage

        da_depth = len(self.obs_attributes)

        np.random.seed(888)
        if da_depth > 1:
            self.obs_noise_seeds = np.random.randint(10000,size=(da_depth)).squeeze()
        else:
            self.obs_noise_seeds = [np.random.randint(10000)]

    @staticmethod
    def converter(results, N, mpv, elem, node, th, ud):
        '''
        Do this after data assimilation for HS balanced vortex.

        '''
        print("Post DA conversion...")

        g = ud.g0
        for n in range(N):
            set_explicit_boundary_data(results[n][0],elem,ud,th,mpv)
            results[n][0].rhoY[...] = (g / 2.0 * results[n][0].rho**2)**th.gamminv

            igy = elem.igy

            kernel = np.ones((2,2))
            kernel /= kernel.sum()

            pn = signal.convolve(results[n][0].rhoY[:,igy,:], kernel, mode='valid')

            set_explicit_boundary_data(results[n][0],elem,ud,th,mpv)
            pn = np.expand_dims(pn, 1)
            pn = np.repeat(pn, node.icy, axis=1)

            results[n][2].p2_nodes[1:-1,:,1:-1] = pn
            set_ghostnodes_p2(results[n][2].p2_nodes,node,ud)

            pn = np.expand_dims(results[n][2].p2_nodes[:,igy,:], 1)
            results[n][2].p2_nodes[...] = np.repeat(pn[...], node.icy, axis=1)

        return results


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
                    data = obs_file[str(attribute)][str(attribute) + str(label)][:]
                    if data.ndim == 3: # implying horizontal slice...
                        data = data[:,0,:]
                    dict_attr = {
                        attribute: data
                    }
                    obs[t_cnt].update(dict_attr)
                t_cnt += 1
            obs = np.array(obs)
            obs_file.close()
        return obs


    def get_loc_mat(self, type, c=1.0):
        x = (np.linspace(0,1,self.obs_X) - 0.5)
        y = (np.linspace(0,1,self.obs_Y) - 0.5)

        X,Y = np.meshgrid(x,y)
        if type == 'constants':
            Z = np.ones((self.obs_X, self.obs_Y))
            return Z.flatten() * c
        
        elif type == 'linear':
            Z = (np.sqrt(X**2 + Y**2))
            Z = Z.max() - Z

            return Z.flatten() * c

        elif type == 'gaussian':
            pos = np.array([X.flatten(),Y.flatten()]).T
            norm = stats.multivariate_normal([0,0],[[0.05,0.0],[0.0,0.05]])
            Z = norm.pdf(pos).reshape(X.shape[0],X.shape[1])
            Z = Z - Z.min()
            Z /= Z.max()

            return Z * c

    def update_dap(self, params):
        for key, value in params.items():
            setattr(self, key, value)


    # setter functions
    @property
    def noise_percentage(self):
        return self._noise_percentage

    @noise_percentage.setter
    def noise_percentage(self, val):
        self._noise_percentage = val
        for key in self.obs_attributes:
            self.obs_noise[key] = self.noise_percentage

    @property
    def obs_attrs(self):
        return self.obs_attributes

    @obs_attrs.setter
    def obs_attrs(self, lst):
        self.obs_attributes = lst
        self.obs_noise = {}
        self.noise_percentage = self._noise_percentage
        if self.add_obs_noise:
            self.gen_obs_noise()
        # setattr(self,self.noise_percentage,self._noise_percentage)

    @property
    def da_times(self):
        return self._da_times

    @da_times.setter
    def da_times(self, lst):
        self._da_times = lst
        self._da_times = np.around(self._da_times,3)
        if self.sparse_obs:
            self.gen_obs_sparse()


    @property
    def loc_setter(self):
        return self.loc_setter

    @loc_setter.setter
    def loc_setter(self, tpl):
        self.obs_X = tpl[0]
        self.obs_Y = tpl[1]
        # constants, linear, gaussian
        self.localisation_matrix = self.get_loc_mat('gaussian')


    @property
    def sd_setter(self):
        return self.sd_setter
    
    @sd_setter.setter
    def sd_setter(self,lst):
        self.std_dev = lst
        self.noise_type = 'FixCov'
        if self.add_obs_noise:
            self.gen_obs_noise()

    # @staticmethod
    # def sampler_gaussian(var):
    #     return lambda ic: np.random.normal(ic,var**0.5)

    # @staticmethod
    # def sampler_none():
    #     return lambda ic: ic

    # @staticmethod
    # def sampler_perturbator(max_shift):
    #     return lambda ic: np.roll(np.roll(ic, np.random.randint(-max_shift,max_shift), axis=1), np.random.randint(-max_shift,max_shift), axis=0)
