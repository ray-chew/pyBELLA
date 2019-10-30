from scipy import sparse, linalg
from scipy.sparse.linalg import spsolve
from scipy.ndimage import map_coordinates
from inputs.enum_bdry import BdryType
import numpy as np

def da_interface(results,obs_current,attr,N,ud,loc=0):
    # inner = (slice(ud.igx,-ud.igx),slice(ud.igy,-ud.igy))
    ig = 2
    inner = (slice(ig,-ig),slice(ig,-ig))
    # inner = (slice(None,),slice(None,))

    local_ens = [getattr(results[:,loc,...][n],attr)[inner] for n in range(N)]
    local_ens = analysis(local_ens,attr)

    # print(obs_current.shape)
    obs_current = obs_current[inner]
    x_obs, y_obs = obs_current.shape
    
    forward_operator = lambda ensemble: interpolation_func(ensemble, x_obs, y_obs, ud)
    local_ens.forward(forward_operator)

    obs_covar = sparse.eye(x_obs**2, y_obs**2, format='csr')
    # print(obs_covar.shape)
    X = local_ens.analyse(obs_current.reshape(-1), obs_covar)
    local_ens.ensemble = local_ens.to_array(X)
    local_ens.ensemble = [np.pad(mem,2,mode='constant') for mem in local_ens.ensemble]
    # print(local_ens.X.shape)
    # return np.array([local_ens.identifier, local_ens.ensemble])
    return [local_ens.identifier, local_ens.ensemble]
    

# let me ust put the forward operator here for now - will need to tidy stuff up....
def interpolation_func(ensemble,x_obs,y_obs,ud):
    if ud.bdry_type[0] == BdryType.WALL or ud.bdry_type[1] == BdryType.WALL:
        assert("WALL NOT IMPLEMENTED!")
    x_ens, y_ens = ensemble[0].shape
    x = np.linspace(0,x_ens,x_obs)
    y = np.linspace(0,y_ens,y_obs)
    x,y = np.meshgrid(x,y)

    ensemble = [map_coordinates(mem,[x,y],mode='wrap', order=3) for mem in ensemble]
    # print(np.array(ensemble).shape)
    return np.array(ensemble)


class analysis(object):
    def __init__(self,ensemble, identifier=None):
        self.ensemble = np.array(ensemble)
        self.X = self.state_vector(ensemble)
        self.no_of_members = self.ensemble.shape[0]
        self.member_shape = self.ensemble[0].shape

        # ensemble inflation factor
        self.rho = 1.3
        # if anlaysis is over local state space, which is it?
        self.identifier = identifier

    def forward(self,forward_operator):
        self.forward_operator = forward_operator

    def localisation(self,localisation_matrix):
        self.localisation_matrix = localisation_matrix

    def analyse(self,obs,obs_covar):
        obs = obs.reshape(-1)
        obs = np.random.normal(obs, obs.max()**0.5)
        #obs: R in l
        #obs_covar: R in (l x l)
        # self.Y = [self.forward_operator(xi) for xi in self.X]
        self.Y = self.forward_operator(self.ensemble)
        self.Y = self.state_vector(self.Y)

        self.Y_mean = self.get_mean(self.Y) # R in l
        self.Y -= self.Y_mean # R in (l x k)

        self.X_mean = self.get_mean(self.X) # R in m
        self.X -= self.X_mean # R in (m x k)

        # for now, global == local, i.e. observation space == state space

        # here is where the R-localisation matrix will come in
        # obs_covar is the error covariance matrix: is spd
        C = spsolve(obs_covar, self.Y.T).T # R in (k x l)

        Pa = (self.no_of_members - 1.) * np.eye(self.no_of_members) / self.rho + np.dot(C,self.Y.T)

        Lambda, P = linalg.eigh(Pa)
        # Lambda, P = Lambda.real, P.real
        # Pa = np.dot(P,np.dot(np.diag(1./Lambda),P.T))
        Pa = P @ np.diag(1./Lambda) @ P.T

        Wa = (self.no_of_members - 1.)**0.5 * np.dot(P,np.dot(np.diag((1./Lambda)**0.5),P.T))

        wa = np.dot(Pa,np.dot(C , (obs - self.Y_mean)))
        # Wa += wa.reshape(1,-1)
        Wa += wa
        # print(Wa.shape)

        return np.dot(self.X.T, Wa.T).T + self.X_mean

    def get_mean(self,vec):
        return np.array(vec).sum(axis=0)
    
    # More readable method needed - seems to be most efficient though.
    @staticmethod
    def state_vector(ensemble):
        return np.array([mem.reshape(-1) for mem in ensemble])
        # return np.array([[getattr(mem,attr).reshape(-1) for mem in ensemble.members(ensemble)] for attr in attributes]).squeeze()

    def to_array(self,X):
        return np.array([x.reshape(self.member_shape) for x in X])