from scipy import sparse, linalg
from scipy.sparse.linalg import spsolve
from scipy.ndimage import map_coordinates
from scipy.signal import convolve2d
from scipy.interpolate import interpn
from inputs.enum_bdry import BdryType
import matplotlib.pyplot as plt
import dask.array as darr
import numpy as np

import data_assimilation.utils as dau

debug_cnt = 0

def da_interface(results,obs_current,rho,attr,N,ud,loc=0):
    # ig = 2
    # inner = (slice(ig,-ig),slice(ig,-ig))
    inner = (slice(None,),slice(None,))

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

    # local_ens.ensemble = np.array([np.pad(mem,ig,mode='constant', constant_values=(0.0)) for mem in local_ens.ensemble])
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


class prepare_rloc(object):
    def __init__(self, ud, elem, node, dap, N, obs_X=5, obs_Y=5):
        # self.elem = elem
        self.iicx = elem.iicx
        self.iicy = elem.iicy
        self.iicxn = node.iicx
        self.iicyn = node.iicy
        self.N = N

        self.i2 = (slice(elem.igx,-elem.igx),slice(elem.igy,-elem.igy))
        self.i0 = (slice(None,),slice(None,))

        self.attr_len = len(dap.obs_attributes)
        self.loc = dap.loc

        self.loc_c = dap.loc_c
        self.loc_f = dap.loc_f
        self.loc_n = dap.loc_n

        self.oa = dap.obs_attributes
        self.da_times = dap.da_times

        self.ca, self.na = self.sort_locs()
        self.cattr_len = len(self.ca)
        self.nattr_len = len(self.na)

        self.obs_X = obs_X
        self.obs_Y = obs_Y

        self.cmask, self.nmask = dau.boundary_mask(ud, elem,node, self.loc_c, self.loc_n)

        self.inf_fac = dap.inflation_factor
        self.loc_mat = dap.localisation_matrix
        

    def sort_locs(self):
        cell_attributes = [key for key,value in self.loc.items() if value == self.loc_c and key in self.oa]
        node_attributes = [key for key,value in self.loc.items() if value == self.loc_n and key in self.oa]
        face_attributes = [key for key,value in self.loc.items() if value == self.loc_f and key in self.oa]

        if len(face_attributes) > 0:
            assert 0, "r-localisation for values on faces not supported."

        return cell_attributes, node_attributes

    def analyse(self,results,obs,N,tout):
        Nx, Ny, obs_attr, attr_len, loc = self.get_properties('cell')

        obs_X, obs_Y = self.obs_X, self.obs_Y

        state_p, obs_p = self.stack(results,obs,obs_attr,tout)

        X = self.get_state(state_p,Nx,Ny,attr_len)
        obs_p = self.get_obs(obs_p,obs_X,obs_Y,Nx,Ny,attr_len)
        Y = self.get_state_in_obs_space(results,obs_attr,obs_X,obs_Y,Nx,Ny,attr_len)

        obs_covar = self.get_obs_covar('cell', Nx, Ny, obs_X, obs_Y)

        analysis_res = np.zeros_like(X)

        for n in range(Nx * Ny):
            # obs_covar_current = sparse.diags(list(obs_covar[n].ravel()) * attr_len, format='csr')
            obs_covar_current = sparse.eye(attr_len*obs_X*obs_Y,attr_len*obs_X*obs_Y, format='csr')

            forward_operator = lambda ensemble : Y[n]
            local_ens = analysis(X[n],self.inf_fac)
            local_ens.forward(forward_operator)
            local_ens.localisation_matrix = self.loc_mat
            analysis_ens = local_ens.analyse(obs_p[n],obs_covar_current)
            analysis_res[n] = analysis_ens

        analysis_res = np.swapaxes(analysis_res,0,2)

        for cnt, attr in enumerate(obs_attr):
            current = analysis_res[cnt]

            for n in range(N):
                data = current[n]
                data = data.reshape(Nx, Ny)
                
                data = np.pad(data,2,mode='constant')

                setattr(results[:,loc,...][n],attr,data)

        return results
    # stack variables on grid according to grid-type. 
    # Currently supports only inner domain with no ghost cells. NO BC handling!
    def stack(self,results,obs,obs_attr,tout):
        state = self.get_quantity(results,obs_attr[0])
        state = state[:,np.newaxis,...]

        obs_stack = np.array(obs[list(self.da_times).index(tout)][obs_attr[0]])[self.i2]
        obs_stack = obs_stack[np.newaxis,...]

        for attr in obs_attr[1:]:
            next_attr = self.get_quantity(results,attr)[:,np.newaxis,...]
            state = np.hstack((state,next_attr))

            next_obs = np.array(obs[list(self.da_times).index(tout)][attr])[self.i2]
            next_obs = next_obs[np.newaxis,...]

            obs_stack = np.vstack((obs_stack,next_obs))

        return state, obs_stack

    def get_state(self,state,Nx,Ny,attr_len):
        X = np.array([dau.sliding_window_view(mem, (1,1), (1,1)).reshape(Nx*Ny,attr_len) for mem in state])
        X = np.swapaxes(X,0,1)
        return X

    # currently supports only 5x5 local counterparts.
    def get_obs(self,obs,obs_X,obs_Y,Nx,Ny,attr_len):
        x_pad = int((obs_X - 1) / 2)
        y_pad = int((obs_Y - 1) / 2)
        obs = np.array([np.pad(obs_attr,(x_pad,y_pad),mode='wrap') for obs_attr in obs])

        obs = np.array([dau.sliding_window_view(obs_attr, (obs_X,obs_Y), (1,1)) for obs_attr in obs])

        obs = np.swapaxes(obs,0,2)
        obs = np.swapaxes(obs,0,1)
        obs = obs.reshape(Nx*Ny,attr_len,obs_X,obs_Y)

        return obs

    def get_state_in_obs_space(self,state,obs_attr,obs_X,obs_Y,Nx,Ny,attr_len):
        sios = self.get_quantity(state,obs_attr[0],inner=False)
        sios = sios[:,np.newaxis,...]

        for attr in obs_attr[1:]:
            next_attr = self.get_quantity(state,attr,inner=False)[:,np.newaxis,...]
            sios = np.hstack((sios,next_attr))

        sios = np.array([dau.sliding_window_view(mem, (obs_X,obs_Y), (1,1)).reshape(Nx*Ny,attr_len,obs_X,obs_Y) for mem in sios])

        sios = np.swapaxes(sios,0,1)
        return sios

    def get_obs_covar(self,type, Nx, Ny, obs_X, obs_Y):
        if type == 'cell':
            obs_covar = np.ones((self.iicx,self.iicy))
        else:
            obs_covar = np.ones((self.iicxn,self.iicyn))

        obs_covar = np.pad(obs_covar, 2, mode='constant', constant_values=(1.0))

        if type == 'cell':
            obs_covar *= self.cmask
        else:
            obs_covar *= self.nmask

        obs_covar = np.array(dau.sliding_window_view(obs_covar, (obs_X,obs_Y), (1,1))).reshape(Nx*Ny,obs_X,obs_Y)

        return obs_covar


    def get_properties(self,type):
        if type != 'cell' and type != 'node':
            assert 0, "rloc: grid-type not supported"

        Nx = self.iicx if type == 'cell' else self.iicxn
        Ny = self.iicy if type == 'cell' else self.iicyn
        obs_attr = self.ca if type == 'cell' else self.na
        attr_len = self.cattr_len if type == 'cell' else self.nattr_len
        loc = self.loc_c if type == 'cell' else self.loc_n

        return Nx, Ny, obs_attr, attr_len, loc
    

    def get_quantity(self,results,quantity,inner=True):
        if inner == True:
            slc = self.i2
        else:
            slc = self.i0

        loc = self.loc[quantity]
        attribute_array = [getattr(results[:,loc,...][n],quantity)[slc] for n in range(self.N)]
        return np.array(attribute_array)



    