import numpy as np

class da_params(object):

    def __init__(self,N):
        # number of ensemble members
        self.N = 20

        # self.attributes = ['rho', 'rhou', 'rhov', 'rhoY']
        self.attributes = ['rho']

        # forward operator (projector from state space to observation space)
        self.forward_operator = np.eye(N)

        # localisation matrix
        self.localisation_matrix = np.eye(N)
        #
        # square of empirical RMSE of (48x48) travelling vortex from ref (256x256)
        self.aprior_error_covar = 0.0001#0.5804227421558537

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