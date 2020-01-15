import numpy as np

class da_params(object):

    def __init__(self,N,da_type='rloc'):
        # number of ensemble members
        self.N = 20
        # self.state_attributes = ['rho', 'rhou', 'rhov','rhow','rhoY','rhoX']
        self.obs_attributes = ['rho', 'rhou', 'rhov']

        # forward operator (projector from state space to observation space)
        self.forward_operator = np.eye(N)

        # localisation matrix
        self.localisation_matrix = np.eye(N)
        weights = [0.05719096,0.25464401,0.33333333,0.25464401,0.05719096,0.25464401,0.52859548,0.66666667,0.52859548,0.25464401,0.33333333,0.66666667,1.,0.66666667,0.33333333,0.25464401,0.52859548,0.66666667,0.52859548,0.25464401,0.05719096,0.25464401,0.33333333,0.25464401,0.05719096]
        # weights = [0.01831564,0.082085,0.13533528,0.082085,0.01831564,0.082085,0.36787944,0.60653066,0.36787944,0.082085,0.13533528,0.60653066,1.,0.60653066,0.13533528,0.082085,0.36787944,0.60653066,0.36787944,0.082085,0.01831564,0.082085,0.13533528,0.082085,0.01831564]
        weights3 = weights*3
        self.localisation_matrix = np.diag(weights3)
        self.localisation_matrix += np.diag(np.ones_like(weights3))
        
        # square of empirical RMSE of (48x48) travelling vortex from ref (256x256)
        self.aprior_error_covar = 0.0001#0.5804227421558537
        self.da_type = da_type

        # ensemble inflation factor for LETKF
        self.inflation_factor = 1.4

        # rejuvenation factor for ETPF
        self.rejuvenation_factor = 0.001

    def load_obs(self,obs):
        None

    @staticmethod
    def sampler_gaussian(var):
        return lambda ic: np.random.normal(ic,var**0.5)

    @staticmethod
    def sampler_none():
        return lambda ic: ic

    @staticmethod
    def sampler_perturbator(max_shift):
        return lambda ic: np.roll(np.roll(ic, np.random.randint(-max_shift,max_shift), axis=1), np.random.randint(-max_shift,max_shift), axis=0)