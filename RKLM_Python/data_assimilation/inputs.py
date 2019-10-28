import numpy as np

class da_params(object):

    def __init__(self,N):
        # number of ensemble members
        self.N = 20

        #
        self.attributes = ['rho']

        # forward operator (projector from state space to observation space)
        self.forward_operator = np.eye(N)

        # localisation matrix
        self.localisation_matrix = np.eye(N)

        #
        self.aprior_error_covar = 0.05

    def load_obs(self,obs):
        None

    @staticmethod
    def sampler_gaussian(var):
        return lambda ic: ic + np.random.normal(0.0,var**0.5,size=ic.shape)

    @staticmethod
    def sampler_none():
        return lambda ic: ic