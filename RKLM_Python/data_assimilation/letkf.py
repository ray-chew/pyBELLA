from scipy import linalg
import numpy as np

class letkf(object):
    def __init__(self,ensemble):
        self.ensemble = ensemble
        self.X = ensemble.state_vector(ensemble)
        self.no_of_members = self.ensemble.members(ensemble).shape[0]

        # ensemble inflation factor
        self.rho = 1.

    def forward(self,forward_operator):
        self.forward_operator = forward_operator.reshape(-1)

    def localisation(self,localisation_matrix):
        self.localisation_matrix = localisation_matrix

    def analyse(self,obs,obs_covar):
        obs = obs.reshape(-1)
        #obs: R in l
        #obs_covar: R in (l x l)
        self.Y = [self.forward_operator * xi for xi in self.X]

        self.Y_mean = self.get_mean(self.Y) # R in l
        self.Y -= self.Y_mean # R in (l x k)

        self.X_mean = self.get_mean(self.X) # R in m
        self.X -= self.X_mean # R in (m x k)

        # for now, global == local, i.e. observation space == state space

        # here is where the R-localisation matrix will come in
        # obs_covar is the error covariance matrix: is spd
        C = linalg.solve(obs_covar, self.Y.T, assume_a='pos').T # R in (k x l)

        Pa = (self.no_of_members - 1.) * np.eye(self.no_of_members) / self.rho + np.dot(C,self.Y.T)

        Lambda, P = linalg.eig(Pa)
        Lambda, P = Lambda.real, P.real
        Pa = P @ np.diag(1./Lambda) @ P.T

        Wa = (self.no_of_members - 1.) * P @ np.diag((1./Lambda)**0.5) @ P.T

        wa = np.dot(np.dot(Pa,C) , obs - self.Y_mean)
        Wa += wa

        return np.dot(self.X.T, Wa).T + self.X_mean

    def get_mean(self,vec):
        mean = np.array(vec).sum(axis=0)
        return mean

    # def get_covariance(self):
    #     self.get_mean()
        
    #     distance = self.ensemble_list - self.ensemble_mean

    #     self.ensemble_covar = 1. / (self.no_of_members - 1.) * (distance.transpose((0,2,1)) @ distance).sum(axis=0)

    #     return self.ensemble_covar

    # def gaussian_salter(self):
    #     # adds Gaussian noise to the ensemble
    #     None
