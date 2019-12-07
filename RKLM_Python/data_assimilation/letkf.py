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

def da_interface(results,obs_current,attr,N,ud,loc=0):
    # inner = (slice(ud.igx,-ud.igx),slice(ud.igy,-ud.igy))
    ig = 2
    inner = (slice(ig,-ig),slice(ig,-ig))
    # inner = (slice(None,),slice(None,))

    local_ens = np.array([getattr(results[:,loc,...][n],attr) for n in range(N)])
    local_ens = np.array([mem[inner] for mem in local_ens])
    print(local_ens.shape)
    local_ens = analysis(local_ens,attr)

    # print(obs_current.shape)
    obs_current = obs_current[inner]
    obs_current = bin_func(obs_current, local_ens.member_shape)
    # local_ens.debug(obs_current,"0before")
    
    x_obs, y_obs = obs_current.shape
    
    # forward_operator = lambda ensemble: interpolation_func(ensemble, x_obs, y_obs, ud)
    forward_operator = lambda ensemble : ensemble
    local_ens.forward(forward_operator)

    obs_covar = sparse.eye(x_obs**2, y_obs**2, format='csr')

    # print(obs_covar.shape)
    X = local_ens.analyse(obs_current.reshape(-1), obs_covar)
    local_ens.ensemble = local_ens.to_array(X)
    # local_ens.debug(obs_current,"1after")
    local_ens.ensemble = np.array([np.pad(mem,2,mode='wrap') for mem in local_ens.ensemble])

    # print(local_ens.X.shape)
    # print(local_ens.ensemble.shape)
    # return np.array([local_ens.identifier, local_ens.ensemble])
    # return [local_ens.identifier, local_ens.ensemble]
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
    # print(np.array(ensemble).shape)
    return np.array(ensemble)


def bin_func(obs,ens_mem_shape):
    obs = obs.reshape(ens_mem_shape[0],obs.shape[0]//ens_mem_shape[0],
                      ens_mem_shape[1],obs.shape[1]//ens_mem_shape[1])
    return obs.mean(axis=(1,3))


def none_func(ensemble):
    return lambda ensemble : ensemble


class analysis(object):
    def __init__(self,ensemble, identifier=None):
        self.ensemble = np.array(ensemble)
        self.X = self.state_vector(ensemble)
        self.no_of_members = self.ensemble.shape[0]
        self.member_shape = self.ensemble[0].shape

        # ensemble inflation factor
        self.rho = 1.4
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
        # obs = np.ones_like(obs) * 100.
        # obs = np.random.normal(obs, np.abs(obs.max()))
        #obs: R in l
        #obs_covar: R in (l x l)
        # self.Y = [self.forward_operator(xi) for xi in self.X]

        self.Y = self.forward_operator(self.X)
        self.Y = self.state_vector(self.Y)

        self.Y_mean = self.get_mean(self.Y) # R in l
        self.Y -= self.Y_mean # R in (l x k)

        self.X_mean = self.get_mean(self.X) # R in m
        self.X -= self.X_mean # R in (m x k)

        # print("self.ensemble.shape = ", self.ensemble.shape)
        # print("self.no_of_members = ", self.no_of_members)
        # print("self.X.shape = ", self.X.shape)
        # print("self.Y.shape = ", self.Y.shape)
        # print("self.Y_mean.shape = ", self.Y_mean.shape)
        # print("obs.shape = ", obs.shape)
        # print("obs_covar.shape = ", obs_covar.shape)
        # print("Sanity check # 1, sum of columns of X = ", self.X.sum(axis=0))
        # col_sum_X = self.X.sum(axis=0)
        # assert(np.allclose(col_sum_X,np.zeros_like(col_sum_X), "sum != 0 "))

        # for now, global == local, i.e. observation space == state space

        # here is where the R-localisation matrix will come in
        # obs_covar is the error covariance matrix: is spd
        # if self.localisation_matrix != None:
            # obs_covar *= self.localisation_matrix
        C = spsolve(obs_covar, self.Y.T).T # R in (k x l)
        if len(self.localisation_matrix) > 0:
            "Starting localisation..."
            # print(self.localisation_matrix.shape)
            # print(C.shape)
            C[...] = ((np.array(self.localisation_matrix)) @ C.T).T
        # print(self.Y.shape)
        # print(obs_covar.shape)
        # C = linalg.solve(obs_covar, self.Y.T).T

        # print("C.shape = ", C.shape)

        # Pa = (self.no_of_members - 1.) * np.eye(self.no_of_members) / self.rho + np.dot(C , self.Y.T)

        Pa = (self.no_of_members - 1.) * np.eye(self.no_of_members) / self.rho + (C @ self.Y.T)

        Lambda, P = linalg.eigh(Pa)
        # Pa = np.dot(P,np.dot(np.diag(1./Lambda),P.T))
        Pa = P @ (np.diag(1./Lambda) @ P.T)

        # print("Pa.shape = ", Pa.shape)

        # Wa = (self.no_of_members - 1.)**0.5 * np.dot(P,np.dot(np.diag((1./Lambda)**0.5),P.T))
        Wa = np.sqrt(self.no_of_members - 1.) * P @ (np.diag(np.sqrt(1./Lambda)) @ P.T)

        # wa = np.dot(Pa , np.dot(C , (obs - self.Y_mean)))
        wa = Pa @ (C @ (obs - self.Y_mean))
        # print(Pa)
        # print(wa)
        # print("Sanity check #2, sum of columns of wa:", np.dot(self.X.T , Wa).sum(axis=1))
        # col_sum = np.dot(self.X.T , Wa).sum(axis=1)
        # assert(np.allclose(col_sum,np.zeros_like(col_sum), "sum != 0 "))

        # Wa += wa.reshape(1,-1)
        Wa += wa
        # print(Wa)

        return (np.dot(self.X.T , Wa.T) + self.X_mean.reshape(-1,1)).T
        # return np.zeros_like(self.X)
        # return ((self.X.T @ Wa) + self.X_mean.reshape(-1,1)).T
        # return result

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