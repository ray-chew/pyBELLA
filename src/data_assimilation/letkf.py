import numpy as np
import numpy.lib.stride_tricks as st

import scipy.sparse as sp
import scipy.ndimage as ndimage

import dask.array as darr
import dask

import matplotlib.pyplot as plt
import logging

import dycore.utils.options as opts
import data_assimilation as da

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

    obs_covar = sp.sparse.eye((x_obs*y_obs), (x_obs*y_obs), format='csr')

    X = local_ens.analyse(obs_current.reshape(-1), obs_covar)
    local_ens.ensemble = local_ens.to_array(X)

    # local_ens.ensemble = np.array([np.pad(mem,ig,mode='constant', constant_values=(0.0)) for mem in local_ens.ensemble])
    # if attr == 'p2_nodes':
    #     local_ens.ensemble = np.array([mem-mem.mean() for mem in local_ens.ensemble])
    # local_ens.ensemble = np.array([np.pad(mem,ig,mode='wrap') for mem in local_ens.ensemble])
    return local_ens.ensemble


# let me ust put the forward operator here for now - will need to tidy stuff up....
def interpolation_func(ensemble,x_obs,y_obs,ud):
    if ud.bdry_type[0] == opts.BdryType.WALL or ud.bdry_type[1] == opts.BdryType.WALL:
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
    ensemble = [ndimage.map_coordinates(mem,[y,x],mode='wrap', order=3) for mem in ensemble]
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

    def __init__(self,ensemble, rho, X_mean, Y_mean, identifier=None):
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

        self.X_mean = X_mean
        self.Y_mean = Y_mean

    def forward(self,forward_operator):
        self.forward_operator = forward_operator

    def localisation(self,localisation_matrix):
        self.localisation_matrix = localisation_matrix

    def analyse(self,obs,obs_covar):
        """
        Analysis method. 'l' is the observation space. 'm' the state space, 'k' the ensemble size.

        """
        if self.forward_operator is None:
            logging.info("Forward operator undefined. Using identity...")
            assert 0, "not implemented."
        if self.localisation_matrix is None:
            logging.info("Localisation matrix undefined. Using identity...")
            # assert(0, "not implemented.")

        # get observations as vector, R in l
        obs = obs.reshape(-1)

        self.Y = self.forward_operator(self.X)
        # get states in observation space as vector, R in (k x l)
        self.Y = self.state_vector(self.Y)

        # self.Y_mean = self.get_mean(self.Y) # R in l
        # self.Y -= self.Y_mean # R in (l x k)

        # self.X_mean = self.get_mean(self.X) # R in m
        # self.X -= self.X_mean # R in (m x k)

        # This is step 4 of the algorithm in Hunt et. al., 2007.
        C = sp.spsolve(obs_covar, self.Y.T).T # R in (k x l)
        # This applies the localisation function to the local region.
        if self.localisation_matrix is not None:
            C[...] = ((np.array(self.localisation_matrix)) @ C.T).T

        # The next three lines are step 5 of the algorithm.
        Pa = (self.no_of_members - 1.) * np.eye(self.no_of_members) / self.rho + (C @ self.Y.T)

        Lambda, P = sp.linalg.eigh(Pa)
        Pa = P @ (np.diag(1./Lambda) @ P.T)

        # This is step 6 of the algorithm.
        Wa = np.sqrt(self.no_of_members - 1.) * P @ (np.diag(np.sqrt(1./Lambda)) @ P.T)

        # The following two lines are step 7 of the algorithm.
        wa = Pa @ (C @ (obs - self.Y_mean))
        Wa += wa

        # This is step 8 of the algorithm.
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

        # get size of the local counterpart.
        self.obs_X = dap.obs_X
        self.obs_Y = dap.obs_Y

        # make sure that subdomain is not larger than inner cellular grid
        assert self.obs_X < elem.iicx, "obs_X > iicx"
        if elem.iicy == 1: # implying horizontal slice
            assert self.obs_Y < elem.iicz, "obs_Y > iicz"
        else:    
            assert self.obs_Y < elem.iicy, "obs_Y > iicy"

        # make sure that subdomain size is odd
        assert self.obs_X % 2 == 1, "obs_X is even"
        assert self.obs_Y % 2 == 1, "obs_Y is even"

        # if checks are passed, then get pad length based on subdomain / local counterpart size
        self.pad_X = int((self.obs_X - 1) / 2)
        self.pad_Y = int((self.obs_Y - 1) / 2)

        # get mask to handle BC. Periodic mask includes ghost cells/nodes and wall masks takes only the inner domain.
        self.cmask, self.nmask = da.utils.boundary_mask(ud, elem, node, self.pad_X, self.pad_Y)

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

    def analyse(self,results,obs,covar,mask,N,tout):
        """
        Wrapper to do analysis by grid-type.

        """

        if self.cattr_len > 0:
            results = self.analyse_by_grid_type(results,obs,covar[0],mask,N,tout,'cell')
        if self.nattr_len > 0:
            results = self.analyse_by_grid_type(results,obs,covar[1],mask,N,tout,'node')
        return results

    def analyse_by_grid_type(self,results,obs,covar,mask,N,tout,grid_type):
        """
        Do analysis by grid-type. Wrapper of the LETKF analysis.

        Note that k is the ensemble size, m is the size of the local state space, attr_len, and l is the size of the local observation space, obs_X*obs*Y*attr_len

        """

        logging.info("Analysis grid type = %s" %grid_type)
        # get properties from class attributes
        Nx, Ny, obs_attr, attr_len, loc = self.get_properties(grid_type)

        # get the size of local counterpart
        obs_X, obs_Y = self.obs_X, self.obs_Y

        # get stacked state space and observation
        # this prepares X, and Y_obs
        state_p, obs_p, mask_p = self.stack(results,obs,mask,obs_attr,tout)

        # get boundary mask
        # this does the BC handling
        bc_mask,mask_n = self.get_bc_mask(mask_p, grid_type, Nx, Ny, obs_X, obs_Y,attr_len)

        # Here, the 2D arrays are split into local counterparts, say of (5x5) arrays. 
        X, X_bar = self.get_state(state_p,mask_p,Nx,Ny,attr_len) # X in R in ((Nx * Ny) x k x m) and X_bar in ((Nx * Ny) * m)
        obs_p = self.get_obs(obs_p,mask_p,obs_X,obs_Y,Nx,Ny,attr_len) # R in ((Nx * Ny) x l)
        Y, Y_bar = self.get_state_in_obs_space(results,mask_p,obs_attr,obs_X,obs_Y,Nx,Ny,attr_len) # Y in R in ((Nx * Ny) x k x l) and Y_bar in ((Nx * Ny) * l)

        analysis_res = np.zeros_like(X)

        # covariance handling
        if covar is None:
            # use identity for observation covariance
            obs_covar_current = sp.sparse.eye(attr_len*obs_X*obs_Y,attr_len*obs_X*obs_Y, format='csr')
        else:
            # get covar at current time
            covar = covar[list(self.da_times).index(tout)]
            
            # expand obs covar to size of the subdomain
            covar = np.expand_dims(covar,axis=-1)
            covar = np.repeat(covar, obs_X, axis=-1)

            covar = np.expand_dims(covar,axis=-1)
            covar = np.repeat(covar, obs_Y, axis=-1)


        # Note: I implemented the easiest and laziest way to make the LETKF work with large localisation regions.
        # However, this link:
        # https://github.com/dask/dask/issues/7589
        # contains a more elegant (and memory efficient) solution using the dask map_overlap function.

        # Calculate chunk size such that largest chunk takes 180mb of memory:
        mem_size = 180
        mb=Y.nbytes/1024/1024 # Y is our largest array
        chunks = int(np.ceil(mb/mem_size))
        chunk_size = int(np.ceil((Nx*Ny) / chunks))

        logging.info("\n===================")
        logging.info("To split DA problem into %i chunks at %s mb each" %(chunks,mem_size))
        logging.info("with analysis of %i grid points per chunk" %chunk_size)
        logging.info("")

        # Now we chunkify our dask arrays:
        bc_mask = darr.rechunk(bc_mask, chunks=chunk_size).to_delayed().ravel()
        mask_n = darr.rechunk(mask_n, chunks=chunk_size).to_delayed().ravel()
        obs_p = darr.rechunk(obs_p, chunks=chunk_size).to_delayed().ravel()
        X = darr.rechunk(X, chunks=chunk_size).to_delayed().ravel()
        X_bar = darr.rechunk(X_bar, chunks=chunk_size).to_delayed().ravel()
        Y = darr.rechunk(Y, chunks=chunk_size).to_delayed().ravel()
        Y_bar = darr.rechunk(Y_bar, chunks=chunk_size).to_delayed().ravel()

        # Now for each chunk, we get the analysis:
        analysis_by_chunk = []
        for chunk in range(chunks):
            analysis_by_chunk.append(darr.from_delayed(self.do_analysis(attr_len, covar, bc_mask[chunk], mask_n[chunk], obs_p[chunk], X[chunk], X_bar[chunk], Y[chunk], Y_bar[chunk]), shape=(attr_len,self.N,chunk_size), dtype=np.int64))

        # Put the results of the chunks together
        analysis_res = darr.concatenate(analysis_by_chunk, axis=2)

        logging.info("Progress of the DA procedure:")
        # Do the computation of the delayed tasks
        with dask.diagnostics.ProgressBar():
            analysis_res = analysis_res.compute()

        logging.info("===================\n")

        # put analysis back into results container
        for cnt, attr in enumerate(obs_attr):
            current = analysis_res[cnt]

            for n in range(N):
                data = current[n]
                data = data.reshape(Nx, Ny)
                
                data = np.pad(data,2,mode='constant')

                setattr(results[:,loc,...][n],attr,data)

        return results


    # loop through all grid-points, selecting either the grid-point or its corresponding local region
    # This is step 3 of the LETKF algorithm in Hunt et. al. 2007. 
    @dask.delayed
    def do_analysis(self, attr_len, covar, bc_mask, mask_n, obs_p, X, X_bar, Y, Y_bar):
        analysis_res = np.zeros_like(X)
        current_chunk_size = analysis_res.shape[0]

        for n in range(current_chunk_size):
            # For each of the quantities in the local observation space, Y, Y_bar, Y_obs, covar and loc_mat, remove the NaNs corresponding to grid-points without observations.

            # using forward operator as a projection of the state into observation space.
            forward_operator = lambda ensemble : np.array([mem[~np.isnan(mem)] for mem in Y[n]])
            Y_mean = Y_bar[n][~np.isnan(Y_bar[n])]

            # setup LETKF class with local state vector
            local_ens = analysis(X[n],self.inf_fac,X_bar[n],Y_mean)

            # setup forward operator method in LETKF class
            local_ens.forward(forward_operator)
            
            # get the localisation matrix with BC handling
            # BC is handled by localisation matrix. Observations on wall ghost cells have zero influence.
            loc_mat = self.get_loc_mat(bc_mask,mask_n,n,attr_len)
            local_ens.localisation_matrix = loc_mat

            # get masked covariance in local domain
            covar_current = np.ma.array(covar,mask=mask_n[n]).filled(fill_value=np.nan)
            covar_current = covar_current[~np.isnan(covar_current)]
            obs_covar_current = sp.sparse.diags(covar_current.flatten(), format='csr')

            # get obs according to sparsity structure
            obs_pn = obs_p[n][~np.isnan(obs_p[n])]

            # do analysis given observations and obs covar.
            analysis_ens = local_ens.analyse(obs_pn,obs_covar_current)

            # This is step 9 of the algorithm, where the grid-points are re-assembled.
            analysis_res[n] = analysis_ens

        analysis_res = np.swapaxes(analysis_res,0,2)
        return analysis_res

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


    def get_state(self,state,mask,Nx,Ny,attr_len):
        """
        Get state vector in grid-localised view. This is step 1 of the LETKF algorithm in Hunt et. al. (2007). Prepares global X and \bar{X}.

        """

        # Get bar{x[G]}, i.e. ensemble averages of the states.
        X_bar = state.mean(axis=0, keepdims=True)

        # Get anomaly for X[G]
        state -= X_bar

        # Decompose X[G] and bar{X[G]} into the local regions view
        # X = np.array([dau.st.sliding_window_view(mem, (1,1), (1,1)).reshape(Nx*Ny,attr_len) for mem in state])
        # X_bar = np.array([dau.st.sliding_window_view(mem, (1,1), (1,1)).reshape(Nx*Ny,attr_len) for mem in X_bar]) 

        state = darr.from_array(state)
        X = st.sliding_window_view(state, (1,1), axis=(2,3)).reshape(self.N,attr_len,Nx*Ny)
        X_bar = darr.from_array(X_bar)
        X_bar = st.sliding_window_view(X_bar, (1,1), axis=(2,3)).reshape(1,attr_len,Nx*Ny)

        X = np.swapaxes(X,1,2)
        X_bar = np.swapaxes(X_bar,1,2)

        # Make first index the local regions
        X = np.swapaxes(X,0,1)
        return X, X_bar[0]

    def get_obs(self,obs,mask,obs_X,obs_Y,Nx,Ny,attr_len):
        """
        Get observations vector in grid-localised view.

        """

        x_pad = int((obs_X - 1) / 2)
        y_pad = int((obs_Y - 1) / 2)
        pads = [[x_pad, x_pad], [y_pad, y_pad]]
        obs = np.array([np.pad(obs_attr,pads,mode='wrap') for obs_attr in obs])
        mask = np.array([np.pad(mask_attr, pads, mode='wrap') for mask_attr in mask])

        obs = np.ma.array(obs,mask=mask).filled(fill_value=np.nan)

        # obs = np.array([dau.st.sliding_window_view(obs_attr, (obs_X,obs_Y), (1,1)) for obs_attr in obs])
        obs = darr.from_array(obs)
        obs = st.sliding_window_view(obs, (obs_X,obs_Y), axis=(1,2))

        obs = np.swapaxes(obs,0,2)
        obs = np.swapaxes(obs,0,1)
        obs = obs.reshape(Nx*Ny,attr_len,obs_X,obs_Y)

        return obs


    def get_state_in_obs_space(self,state,mask,obs_attr,obs_X,obs_Y,Nx,Ny,attr_len):
        """
        Get state vector projected onto observation space, in grid-localised view. This is step 1 of the LETKF algorithm in Hunt et. al., 2007.

        """

        # Prepare data structure for Y[G] calculation
        # stack grid by attributes 
        sios = self.get_quantity(state,obs_attr[0],inner=False)

        sios = sios[:,np.newaxis,...]

        for attr in obs_attr[1:]:
            next_attr = self.get_quantity(state,attr,inner=False)
            # next_attr = np.ma.array(next_attr,mask=mask).filled(fill_value=np.nan)
            next_attr = next_attr[:,np.newaxis,...]
            sios = np.hstack((sios,next_attr))

        # prepare mask array
        # grid-points where there are observations is 0
        x_pad = int((obs_X - 1) / 2)
        y_pad = int((obs_Y - 1) / 2)
        pads = [[x_pad, x_pad], [y_pad, y_pad]]
        mask = np.array([np.pad(mask_attr, pads, mode='wrap') for mask_attr in mask])

        # apply mask to get Y[G]. Now Y[G] is in observation space. 
        sios = np.ma.array([np.ma.array(mem,mask=mask) for mem in sios])

        # Get bar{y[G]}, i.e. ensemble averages of the states in obs space.
        Y_bar = sios.mean(axis=0,keepdims=True)
        # Get anomaly for Y[G]
        sios -= Y_bar

        # Replace masked values by NaNs for decomposition
        sios = sios.filled(np.nan)
        Y_bar = Y_bar.filled(np.nan)

        # Decompose Y[G] and bar{y[G]} such that the local regions can be accessed easily
        # sios = np.array([st.sliding_window_view(mem, (obs_X,obs_Y), axis=(1,2)).reshape(Nx*Ny,attr_len,obs_X,obs_Y) for mem in sios])
        # Y_bar = np.array([st.sliding_window_view(mem, (obs_X,obs_Y), axis=(1,2)).reshape(Nx*Ny,attr_len,obs_X,obs_Y) for mem in Y_bar])

        sios = darr.from_array(sios)
        sios = st.sliding_window_view(sios, (obs_X,obs_Y), axis=(2,3)).reshape(self.N,attr_len,Nx*Ny,obs_X,obs_Y)
        Y_bar = darr.from_array(Y_bar)
        Y_bar = st.sliding_window_view(Y_bar, (obs_X,obs_Y), axis=(2,3)).reshape(1,attr_len,Nx*Ny,obs_X,obs_Y)

        ###############################

        # from skimage.util import view_as_windows
        # sios_tmp = np.zeros((10,attr_len,Nx*Ny,obs_X,obs_Y))
        # Y_bar_tmp = np.zeros((1,attr_len,Nx*Ny,obs_X,obs_Y))
        
        # # sios_tmp = darr.from_array(sios_tmp)
        # sios = darr.from_array(sios)

        # hN = int(np.floor(self.N / 2))
        # # print(hN)
        # tmp0 = deepcopy(sios[:hN,...])
        # tmp1 = deepcopy(sios[hN:self.N,...])

        # for idx in range(0,self.N):
        #     idx = (self.N - 1) - idx
        #     print(idx)
        #     sios_tmp[idx] = st.sliding_window_view(sios[idx], (obs_X,obs_Y),axis=(1,2)).reshape(attr_len,Nx*Ny,obs_X,obs_Y)


        # sios_tmp00 = np.zeros((5,attr_len,Nx*Ny,obs_X,obs_Y))

        # for idx in range(0,self.N-hN):
        #     print(idx+hN)
        #     sios_tmp00[idx] = st.sliding_window_view(sios[idx+hN], (obs_X,obs_Y),axis=(1,2)).reshape(attr_len,Nx*Ny,obs_X,obs_Y)
        

        ####################################

        #     # sios_tmp[idx] = view_as_windows(mem, (attr_len,obs_X,obs_Y)).reshape(attr_len,Nx*Ny,obs_X,obs_Y)
        # Y_bar_tmp[0] = view_as_windows(Y_bar[0], (attr_len,obs_X,obs_Y)).reshape(attr_len,Nx*Ny,obs_X,obs_Y)
        # sios_tmp = view_as_windows(sios, (self.N,attr_len,obs_X,obs_Y)).reshape(self.N,Nx*Ny,attr_len,obs_X,obs_Y)

        # Y_bar_tmp[0] = st.sliding_window_view(Y_bar[0], (obs_X,obs_Y),axis=(1,2)).reshape(attr_len,Nx*Ny,obs_X,obs_Y)
        
        
        # sios = np.array([st.sliding_window_view(mem, (obs_X,obs_Y), axis=(1,2)).reshape(attr_len,Nx*Ny,obs_X,obs_Y) for mem in sios])
        # Y_bar = np.array([st.sliding_window_view(mem, (obs_X,obs_Y), axis=(1,2)).reshape(attr_len,Nx*Ny,obs_X,obs_Y) for mem in Y_bar])


        # sios = np.array([view_as_windows(mem, (attr_len, obs_X,obs_Y)).reshape(attr_len,Nx*Ny,obs_X,obs_Y) for mem in sios])
        # Y_bar = np.array([view_as_windows(mem, (attr_len, obs_X,obs_Y)).reshape(attr_len,Nx*Ny,obs_X,obs_Y) for mem in Y_bar])

        # sios = sios_tmp
        # Y_bar = Y_bar_tmp


        # Y_bar = np.array([dau.st.sliding_window_view(mem, (obs_X,obs_Y), (1,1)).reshape(Nx*Ny,attr_len,obs_X,obs_Y) for mem in Y_bar])

        sios = np.swapaxes(sios,1,2)
        Y_bar = np.swapaxes(Y_bar,1,2)

        # Make local region the first axis index
        sios = np.swapaxes(sios,0,1)

        # sios = sios.rechunk(chunks=(self.N,100,attr_len,obs_X,obs_Y))
        # print(sios.chunks)
        return sios, Y_bar[0]


    def get_bc_mask(self, mask, type, Nx, Ny, obs_X, obs_Y,attr_len):
        """
        Mask handling the boundary condition for cell or node grids in the local subdomains.

        """
        if type == 'cell' and self.iicy > 1:
            bc_mask = np.ones((self.iicx,self.iicy))
        elif type == 'cell' and self.iicy == 1:
            bc_mask = np.ones((self.iicx,self.iicz))
        elif type == 'node' and self.iicyn > 2:
            bc_mask = np.ones((self.iicxn,self.iicyn))
        elif type == 'node' and self.iicyn == 2:
            bc_mask = np.ones((self.iicxn,self.iiczn))

        # bc_mask = np.pad(bc_mask, 2, mode='constant', constant_values=(1.0))

        pads = ((self.pad_X,self.pad_X),(self.pad_Y,self.pad_Y))

        bc_mask = np.pad(bc_mask, pads, mode='constant', constant_values=(1.0))

        mask = np.array([np.pad(attr, pads, mode='wrap') for attr in mask])

        if type == 'cell':
            bc_mask *= self.cmask
        else:
            bc_mask *= self.nmask

        # bc_mask = np.array(dau.st.sliding_window_view(bc_mask, (obs_X,obs_Y), (1,1))).reshape(Nx*Ny,obs_X,obs_Y)

        # mask_n = np.array(dau.st.sliding_window_view(mask, (obs_X,obs_Y), (1,1))).reshape(Nx*Ny,attr_len,obs_X,obs_Y)

        bc_mask = darr.from_array(bc_mask)
        bc_mask = st.sliding_window_view(bc_mask, (obs_X,obs_Y)).reshape(Nx*Ny,obs_X,obs_Y)

        mask = darr.from_array(mask)
        mask_n = st.sliding_window_view(mask, (obs_X,obs_Y), axis=(1,2)).reshape(attr_len,Nx*Ny,obs_X,obs_Y)

        mask_n = np.swapaxes(mask_n,0,1)

        return bc_mask, mask_n


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
            assert 0, "rloc: grid-type not supported"

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
        slc = self.i2
        loc = self.loc[quantity]
        pads = ((self.pad_X,self.pad_X),(self.pad_Y,self.pad_Y))

        if inner:
            attribute_array = [getattr(results[:,loc,...][n],quantity)[slc] for n in range(self.N)]
        else:
            attribute_array = [np.pad(getattr(results[:,loc,...][n],quantity)[slc],pads, mode='wrap') for n in range(self.N)]

        return np.array(attribute_array)


    def get_loc_mat(self,bc_mask,mask_n,n,attr_len):
        loc_mat = self.loc_mat * bc_mask[n]

        loc_mat = np.expand_dims(loc_mat,axis=0)
        loc_mat = np.repeat(loc_mat, attr_len, axis=0)

        loc_mat = np.ma.array(loc_mat,mask=mask_n[n]).filled(fill_value=np.nan)

        loc_mat = loc_mat[~np.isnan(loc_mat)]
        loc_mat = np.diag(loc_mat.flatten())
        return loc_mat


    
