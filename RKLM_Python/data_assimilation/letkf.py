from scipy import sparse, linalg
from scipy.sparse.linalg import spsolve
from scipy.ndimage import map_coordinates
from scipy.signal import convolve2d
from scipy.interpolate import interpn
from inputs.enum_bdry import BdryType
import matplotlib.pyplot as plt
import dask.array as darr
import numpy as np

debug_cnt = 0

def da_interface(results,obs_current,rho,attr,N,ud,loc=0):
    ig = 2
    inner = (slice(ig,-ig),slice(ig,-ig))
    # inner = (slice(None,),slice(None,))

    local_ens = np.array([getattr(results[:,loc,...][n],attr) for n in range(N)])
    local_ens = np.array([mem[inner] for mem in local_ens])
    local_ens = analysis(local_ens,rho,attr)

    obs_current = obs_current[inner]
    # obs_current = bin_func(obs_current, local_ens.member_shape)
    # local_ens.debug(obs_current,"0before")
    
    x_obs, y_obs = obs_current.shape
    
    forward_operator = lambda ensemble : ensemble
    local_ens.forward(forward_operator)

    obs_covar = sparse.eye((x_obs*y_obs), (x_obs*y_obs), format='csr')

    X = local_ens.analyse(obs_current.reshape(-1), obs_covar)
    local_ens.ensemble = local_ens.to_array(X)

    local_ens.ensemble = np.array([np.pad(mem,ig,mode='constant', constant_values=(0.0)) for mem in local_ens.ensemble])
    # if attr == 'p2_nodes':
    #     local_ens.ensemble = np.array([mem-mem.mean() for mem in local_ens.ensemble])
    # local_ens.ensemble = np.array([np.pad(mem,ig,mode='wrap') for mem in local_ens.ensemble])
    return local_ens.ensemble


# let me ust put the forward operator here for now - will need to tidy stuff up....
def interpolation_func(ensemble,x_obs,y_obs,ud):
    if ud.bdry_type[0] == BdryType.WALL or ud.bdry_type[1] == BdryType.WALL:
        assert("WALL NOT IMPLEMENTED!")
    x_ens, y_ens = ensemble[0].shape
    x = np.linspace(0,x_ens,x_obs)
    y = np.linspace(0,y_ens,y_obs)

    # x = np.linspace(-0.5,0.5,x_ens)
    # y = np.linspace(-0.5,0.5,y_ens)

    # x0 = np.linspace(-0.5,0.5,x_obs)
    # y0 = np.linspace(-0.5,0.5,y_obs)
    # mesh = np.array(np.meshgrid(x0,y0))n
    # pts = np.rollaxis(mesh, 0, 3).reshape((-1, 2))

    # ensemble = [interpn((x,y),mem,pts, method='splinef2d').reshape(x_obs,y_obs) for mem in ensemble]

    x,y = np.meshgrid(x,y)
    ensemble = [map_coordinates(mem,[y,x],mode='wrap', order=3) for mem in ensemble]

    return np.array(ensemble)


def bin_func(obs,ens_mem_shape):
    obs = obs.reshape(ens_mem_shape[0],obs.shape[0]//ens_mem_shape[0],
                      ens_mem_shape[1],obs.shape[1]//ens_mem_shape[1])
    return obs.mean(axis=(1,3))


def none_func(ensemble):
    return lambda ensemble : ensemble


class analysis(object):
    def __init__(self,ensemble, rho, identifier=None):
        self.ensemble = np.array(ensemble)
        self.X = self.state_vector(ensemble)
        self.no_of_members = self.ensemble.shape[0]
        self.member_shape = self.ensemble[0].shape

        # ensemble inflation factor
        self.rho = rho
        # if anlaysis is over local state space, which is it?
        self.identifier = identifier

        # initialise localisation matrix as None
        self.localisation_matrix = []

    def forward(self,forward_operator):
        self.forward_operator = forward_operator

    def localisation(self,localisation_matrix):
        self.localisation_matrix = localisation_matrix

    def analyse(self,obs,obs_covar):
        obs = obs.reshape(-1)

        self.Y = self.forward_operator(self.X)
        self.Y = self.state_vector(self.Y)

        self.Y_mean = self.get_mean(self.Y) # R in l
        self.Y -= self.Y_mean # R in (l x k)

        self.X_mean = self.get_mean(self.X) # R in m
        self.X -= self.X_mean # R in (m x k)

        C = spsolve(obs_covar, self.Y.T).T # R in (k x l)
        if len(self.localisation_matrix) > 0:
            C[...] = ((np.array(self.localisation_matrix)) @ C.T).T

        Pa = (self.no_of_members - 1.) * np.eye(self.no_of_members) / self.rho + (C @ self.Y.T)

        Lambda, P = linalg.eigh(Pa)
        Pa = P @ (np.diag(1./Lambda) @ P.T)

        Wa = np.sqrt(self.no_of_members - 1.) * P @ (np.diag(np.sqrt(1./Lambda)) @ P.T)
        wa = Pa @ (C @ (obs - self.Y_mean))

        Wa += wa

        return (np.dot(self.X.T , Wa.T) + self.X_mean.reshape(-1,1)).T


    def get_mean(self,vec):
        return np.array(vec).mean(axis=0)
    
    # More readable method needed - seems to be most efficient though.
    @staticmethod
    def state_vector(ensemble):
        return np.array([mem.reshape(-1) for mem in ensemble])
        # return np.array([[getattr(mem,attr).reshape(-1) for mem in ensemble.members(ensemble)] for attr in attributes]).squeeze()

    def to_array(self,X):
        return np.array([x.reshape(self.member_shape) for x in X])

    def debug(self, obs_current, suffix=""):
        global debug_cnt
        N = self.ensemble.shape[0]
        ens_mean = self.ensemble.mean(axis=0)
        _, ax = plt.subplots(ncols=N+1, figsize=(15,4))
        for n in range(N):
            ax[n].imshow(self.ensemble[n])
        ax[n+1].imshow(obs_current)
        plt.savefig("./output_images/da_debug_%s_%i_%s" %(self.identifier,debug_cnt // 3,suffix), bbox_inches='tight')
        plt.close()
        debug_cnt += 1