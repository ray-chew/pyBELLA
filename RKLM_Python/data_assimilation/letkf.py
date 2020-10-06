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

def da_interface(results,dap,obs,attr,tout,N,ud):
    """
    Interface for batch-observations localisation with LETKF.

    """

    obs_current = np.array(obs[list(dap.da_times).index(tout)][attr])
    inf_fac = dap.inflation_factor
    loc = dap.loc[attr]

    inner = (slice(None,),slice(None,))

    local_ens = np.array([getattr(results[:,loc,...][n],attr) for n in range(N)])
    local_ens = np.array([mem[inner] for mem in local_ens])
    local_ens = analysis(local_ens,inf_fac,attr)

    obs_current = obs_current[inner]
    
    x_obs, y_obs = obs_current.shape
    
    forward_operator = lambda ensemble : ensemble # state space == observation space
    # forward_operator = dap.forward_operator
    # localisation_matrix = dap.localisation_matrix
    # localisation_matrix = np.ones_like(obs_current.flatten())
    local_ens.forward(forward_operator)
    # local_ens.localisation(localisation_matrix)

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
    """
    LETKF Analysis based on (Hunt et al., 2007). 

    """

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
        self.forward_operator = None
        self.localisation_matrix = None

    def forward(self,forward_operator):
        self.forward_operator = forward_operator

    def localisation(self,localisation_matrix):
        self.localisation_matrix = localisation_matrix

    def analyse(self,obs,obs_covar):
        """
        Analysis method. 'l' is the observation space. 'm' the state space, 'k' the ensemble size.

        """
        if self.forward_operator is None:
            print("Forward operator undefined. Using identity...")
            assert(0, "not implemented.")
        if self.localisation_matrix is None:
            print("Localisation matrix undefined. Using identity...")
            # assert(0, "not implemented.")

        obs = obs.reshape(-1)

        self.Y = self.forward_operator(self.X)
        self.Y = self.state_vector(self.Y)

        self.Y_mean = self.get_mean(self.Y) # R in l
        self.Y -= self.Y_mean # R in (l x k)

        self.X_mean = self.get_mean(self.X) # R in m
        self.X -= self.X_mean # R in (m x k)

        C = spsolve(obs_covar, self.Y.T).T # R in (k x l)
        if self.localisation_matrix is not None:
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
    """
    Helper class to get grid-based localisation for LETKF. Used only when da_type=='rloc' is True.
    
    """
    def __init__(self, ud, elem, node, dap, N, obs_X=5, obs_Y=5):
        # get grid properties
        self.iicx = elem.iicx
        self.iicy = elem.iicy
        self.iicz = elem.iicz
        self.iicxn = node.iicx
        self.iicyn = node.iicy
        self.iiczn = node.iicz
        self.N = N

        # define inner and outer domains
        self.i2 = (slice(elem.igx,-elem.igx),slice(elem.igy,-elem.igy))
        self.i0 = (slice(None,),slice(None,))

        # get da parameters
        self.attr_len = len(dap.obs_attributes)
        self.loc = dap.loc

        self.loc_c = dap.loc_c
        self.loc_f = dap.loc_f
        self.loc_n = dap.loc_n

        self.oa = dap.obs_attributes
        self.da_times = dap.da_times

        # get cell and node-based attributes for DA
        self.ca, self.na = self.sort_locs()
        self.cattr_len = len(self.ca)
        self.nattr_len = len(self.na)

        # get size of the local counterpart. default = (5x5).
        self.obs_X = obs_X
        self.obs_Y = obs_Y

        # get mask to handle BC. Periodic mask includes ghost cells/nodes and wall masks takes only the inner domain.
        self.cmask, self.nmask = dau.boundary_mask(ud, elem, node)

        # get from da parameters the localisation matrix and inflation factor
        self.inf_fac = dap.inflation_factor
        self.loc_mat = dap.localisation_matrix
        

    def sort_locs(self):
        """
        Given a dictionary of attributes and the index (locs) of its data container, return attributes that are cell-based and node-based.

        Note
        ----
        Face-based flux is not yet supported.

        """

        cell_attributes = [key for key,value in self.loc.items() if value == self.loc_c and key in self.oa]
        node_attributes = [key for key,value in self.loc.items() if value == self.loc_n and key in self.oa]
        face_attributes = [key for key,value in self.loc.items() if value == self.loc_f and key in self.oa]

        if len(face_attributes) > 0:
            assert 0, "r-localisation for values on faces not supported."

        return cell_attributes, node_attributes

    def analyse(self,results,obs,mask,N,tout):
        """
        Wrapper to do analysis by grid-type.

        """

        if self.cattr_len > 0:
            results = self.analyse_by_grid_type(results,obs,mask,N,tout,'cell')
        if self.nattr_len > 0:
            results = self.analyse_by_grid_type(results,obs,mask,N,tout,'node')
        return results

    def analyse_by_grid_type(self,results,obs,mask,N,tout,grid_type):
        """
        Do analysis by grid-type. Wrapper of the LETKF analysis.

        """

        print("Analysis grid type = %s" %grid_type)
        # get properties from class attributes
        Nx, Ny, obs_attr, attr_len, loc = self.get_properties(grid_type)

        # get the size of local counterpart
        obs_X, obs_Y = self.obs_X, self.obs_Y

        # get stacked state space and observation
        state_p, obs_p, mask_p = self.stack(results,obs,mask,obs_attr,tout)

        mask_p = mask_p[0].flatten()

        # Here, the 2D arrays are split into local counterparts, say of (5x5) arrays. 
        X = self.get_state(state_p,Nx,Ny,attr_len)
        obs_p = self.get_obs(obs_p,obs_X,obs_Y,Nx,Ny,attr_len)
        Y = self.get_state_in_obs_space(results,obs_attr,obs_X,obs_Y,Nx,Ny,attr_len)

        bc_mask = self.get_bc_mask(grid_type, Nx, Ny, obs_X, obs_Y)

        analysis_res = np.zeros_like(X)

        # loop through all grid-points
        for n in range(Nx * Ny):
            if mask_p[n]:
                # if no observation is at grid location, analysis = forecast.
                analysis_res[n] = X[n]
            else:
                # using identity for observation covariance
                obs_covar_current = sparse.eye(attr_len*obs_X*obs_Y,attr_len*obs_X*obs_Y, format='csr')

                # using forward operator as a projection of the state into observation space.
                forward_operator = lambda ensemble : Y[n]

                # setup LETKF class with state vector
                local_ens = analysis(X[n],self.inf_fac)

                # setup forward operator method in LETKF class
                local_ens.forward(forward_operator)
                
                # get the localisation matrix with BC handling
                # BC is handled by localisation matrix.
                local_ens.localisation_matrix = np.diag(list(bc_mask[n].ravel()) * attr_len) * np.diag((list(self.loc_mat) * attr_len))

                # do analysis given observations and obs covar.
                analysis_ens = local_ens.analyse(obs_p[n],obs_covar_current)

                analysis_res[n] = analysis_ens

        analysis_res = np.swapaxes(analysis_res,0,2)

        # put analysis back into results container
        for cnt, attr in enumerate(obs_attr):
            current = analysis_res[cnt]

            for n in range(N):
                data = current[n]
                data = data.reshape(Nx, Ny)
                
                data = np.pad(data,2,mode='constant')

                setattr(results[:,loc,...][n],attr,data)

        return results


    def stack(self,results,obs,mask,obs_attr,tout):
        """
        On each grid-point, stack all the quantities to be assimilated. This stacking is done separately for cells and nodes.

        Quantities, e.g. {rho, rhou, ...}.

        """
        state = self.get_quantity(results,obs_attr[0])
        state = state[:,np.newaxis,...]

        obs_stack = np.array(obs[list(self.da_times).index(tout)][obs_attr[0]])[self.i2]
        obs_stack = obs_stack[np.newaxis,...]

        mask_stack = np.array(mask[list(self.da_times).index(tout)][obs_attr[0]])[self.i2]
        mask_stack = mask_stack[np.newaxis,...]

        for attr in obs_attr[1:]:
            next_attr = self.get_quantity(results,attr)[:,np.newaxis,...]
            state = np.hstack((state,next_attr))

            next_obs = np.array(obs[list(self.da_times).index(tout)][attr])[self.i2]
            next_obs = next_obs[np.newaxis,...]

            obs_stack = np.vstack((obs_stack,next_obs))

            next_mask = np.array(mask[list(self.da_times).index(tout)][attr])[self.i2]
            next_mask = next_mask[np.newaxis,...]

            mask_stack = np.vstack((mask_stack,next_mask))

        return state, obs_stack, mask_stack


    def get_state(self,state,Nx,Ny,attr_len):
        """
        Get state vector in grid-localised view.

        """

        X = np.array([dau.sliding_window_view(mem, (1,1), (1,1)).reshape(Nx*Ny,attr_len) for mem in state])
        X = np.swapaxes(X,0,1)
        return X


    def get_obs(self,obs,obs_X,obs_Y,Nx,Ny,attr_len):
        """
        Get observations vector in grid-localised view.

        """

        x_pad = int((obs_X - 1) / 2)
        y_pad = int((obs_Y - 1) / 2)
        obs = np.array([np.pad(obs_attr,(x_pad,y_pad),mode='wrap') for obs_attr in obs])

        obs = np.array([dau.sliding_window_view(obs_attr, (obs_X,obs_Y), (1,1)) for obs_attr in obs])

        obs = np.swapaxes(obs,0,2)
        obs = np.swapaxes(obs,0,1)
        obs = obs.reshape(Nx*Ny,attr_len,obs_X,obs_Y)

        return obs


    def get_state_in_obs_space(self,state,obs_attr,obs_X,obs_Y,Nx,Ny,attr_len):
        """
        Get state vector projected onto observation space, in grid-localised view.

        """
        sios = self.get_quantity(state,obs_attr[0],inner=False)
        sios = sios[:,np.newaxis,...]

        for attr in obs_attr[1:]:
            next_attr = self.get_quantity(state,attr,inner=False)[:,np.newaxis,...]
            sios = np.hstack((sios,next_attr))

        sios = np.array([dau.sliding_window_view(mem, (obs_X,obs_Y), (1,1)).reshape(Nx*Ny,attr_len,obs_X,obs_Y) for mem in sios])

        sios = np.swapaxes(sios,0,1)
        return sios


    def get_bc_mask(self,type, Nx, Ny, obs_X, obs_Y):
        """
        Mask handling the boundary condition for either cell or node grids.

        """
        if type == 'cell' and self.iicy > 1:
            bc_mask = np.ones((self.iicx,self.iicy))
        elif type == 'cell' and self.iicy == 1:
            bc_mask = np.ones((self.iicx,self.iicz))
        elif type == 'node' and self.iicyn > 2:
            bc_mask = np.ones((self.iicxn,self.iicyn))
        elif type == 'node' and self.iicyn == 2:
            bc_mask = np.ones((self.iicxn,self.iiczn))

        bc_mask = np.pad(bc_mask, 2, mode='constant', constant_values=(1.0))

        if type == 'cell':
            bc_mask *= self.cmask
        else:
            bc_mask *= self.nmask

        bc_mask = np.array(dau.sliding_window_view(bc_mask, (obs_X,obs_Y), (1,1))).reshape(Nx*Ny,obs_X,obs_Y)

        return bc_mask


    def get_properties(self,type):
        """
        For a given grid-type (cell / node), return the 2D grid-size (Nx,Ny), the attributes {rho, rhou...} observed on this grid (obs_attr), the number of attributes (attr_len), and the index location of its data container.

        Returns
        -------
        Nx : int
        Ny : int
        obs_attr : list of str
        attr_len : int
        loc : int

        """
        if type != 'cell' and type != 'node':
            assert(0, "rloc: grid-type not supported")

        Nx = self.iicx if type == 'cell' else self.iicxn
        if self.iicy > 1:
            Ny = self.iicy if type == 'cell' else self.iicyn
        else:
            Ny = self.iicz if type == 'cell' else self.iiczn
        obs_attr = self.ca if type == 'cell' else self.na
        attr_len = self.cattr_len if type == 'cell' else self.nattr_len
        loc = self.loc_c if type == 'cell' else self.loc_n

        return Nx, Ny, obs_attr, attr_len, loc
    

    def get_quantity(self,results,quantity,inner=True):
        """
        Get ensemble representation of {rho, rhou...}.

        """
        if inner == True:
            slc = self.i2
        else:
            slc = self.i0

        loc = self.loc[quantity]
        attribute_array = [getattr(results[:,loc,...][n],quantity)[slc] for n in range(self.N)]
        return np.array(attribute_array)



    