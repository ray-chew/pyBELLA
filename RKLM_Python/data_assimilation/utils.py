import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from scipy import signal
from inputs import boundary, enum_bdry

class ensemble(object):
    def __init__(self, input_ensemble=[None]):
        if input_ensemble[0] != None:
            cnt = 0
            for mem in input_ensemble:
                setattr(self,'mem_' + str(cnt),mem)
                # self.debug_im(mem[0].rho, cnt)
                cnt += 1
        else:
            None
        

    def initialise_members(self,ic,N):
        for cnt in range(N):
            mem = [deepcopy(arr) for arr in ic]
            # mem = sampler(mem)
            setattr(self,'mem_' + str(cnt),mem)

    # def state_vector(self,ensemble):
    #     N = self.members(ensemble).shape[0]
    #     return self.members(ensemble).reshape(N,-1)

    def set_members(self,analysis_ensemble, tout):
        cnt = 0
        for xi in analysis_ensemble:
            setattr(self,'mem_' + str(cnt),np.array(xi))
            self.debug_im(xi[0].rho, cnt, tout)
            cnt += 1

    # rethink this eveutally....
    def ensemble_spreading(self, ens, sampler, attributes, loc=0):
        N = self.members(ens).shape[0]
        np.random.seed(555)
        for attribute in attributes:
            for n in range(N):
                mem = getattr(self,'mem_' + str(n))
                value = getattr(mem[loc],attribute)
                # self.debug_im(sampler(value), n)
                setattr(mem[loc],attribute,sampler(value))

    @staticmethod
    def members(ensemble):
        return np.array(list(ensemble.__dict__.values()))

    @staticmethod
    def debug_im(value, n, tout):
        plt.figure()
        plt.imshow(value, origin='lower')
        plt.savefig("./output_images/ensemble_snapshots/%.3f_%03d.png" %(tout,n), bbox_inches='tight')
        plt.close()

def ensemble_inflation(results, attributes, factor, N, loc=0):
    for attribute in attributes:
        mean = [getattr(results[n][loc],attribute) for n in range(N)]
        mean = np.array(mean)

        mean = np.mean(mean,axis=0)
        for n in range(N):
            inflation = mean + factor * (getattr(results[n][loc],attribute) - mean)
            setattr(results[n][loc],attribute,inflation)

def set_p2_nodes(analysis,results,N,th,node,ud,loc_c=0,loc_n=2):
    for n in range(N):
        rhoY = analysis[n]
        p2_n = getattr(results[n][loc_n],'p2_nodes')
        rhoY_n = np.zeros_like(p2_n)
        kernel = np.array([[1.,1.],[1.,1.]])
        rhoY_n[1:-1,1:-1] = signal.fftconvolve(rhoY,kernel,mode='valid') / kernel.sum()
        boundary.set_ghostnodes_p2(rhoY_n,node,ud)
        p2_n = rhoY_n**th.gm1 - 1.0 + (p2_n - p2_n.mean())
        # p2_n = rhoY_n**th.gm1 - 1.0
        p2_n -= p2_n.mean()
        # p2_n = np.pad(p2_n,2,mode='wrap')
        boundary.set_ghostnodes_p2(p2_n,node,ud)
        setattr(results[n][loc_n],'p2_nodes',p2_n)

def set_rhoY_cells(analysis,results,N,th,ud,loc_c=0,loc_n=2):
    for n in range(N):
        p2n = analysis[n]
        rhoYc0 = getattr(results[n][loc_c], 'rhoY')
        kernel = np.array([[1.,1.],[1.,1.]])
        p2c = signal.fftconvolve(p2n,kernel,mode='valid') / kernel.sum()
        p2c -= p2c.mean()

        rhoYc = (rhoYc0**th.gm1 + ud.Msq*p2c)**th.gm1inv
        setattr(results[n][loc_c], 'rhoYc', rhoYc)
        
def boundary_mask(ud,elem,node,loc_c,loc_n):
    cmask = np.ones(elem.isc).squeeze()
    nmask = np.ones(node.isc).squeeze()

    for dim in range(elem.ndim):
        ghost_padding = [[0,0]] * elem.ndim
        ghost_padding[dim] = [elem.igs[dim],elem.igs[dim]]

        if ud.bdry_type[dim] == enum_bdry.BdryType.PERIODIC:
            cmask = np.pad(cmask, ghost_padding, mode='constant', constant_values=(1.0))
            nmask = np.pad(nmask, ghost_padding, mode='constant', constant_values=(1.0))

        elif ud.bdry_type[dim] == enum_bdry.BdryType.WALL:
            cmask = np.pad(cmask, ghost_padding, mode='constant', constant_values=(0.0))
            nmask = np.pad(nmask, ghost_padding, mode='constant', constant_values=(0.0))
    
    return cmask, nmask





# ref: https://gist.github.com/meowklaski/4bda7c86c6168f3557657d5fb0b5395a
def sliding_window_view(arr, window_shape, steps):
    """ Produce a view from a sliding, striding window over `arr`.
        The window is only placed in 'valid' positions - no overlapping
        over the boundary.
        Parameters
        ----------
        arr : numpy.ndarray, shape=(...,[x, (...), z])
            The array to slide the window over.
        window_shape : Sequence[int]
            The shape of the window to raster: [Wx, (...), Wz],
            determines the length of [x, (...), z]
        steps : Sequence[int]
            The step size used when applying the window
            along the [x, (...), z] directions: [Sx, (...), Sz]
        Returns
        -------
        view of `arr`, shape=([X, (...), Z], ..., [Wx, (...), Wz])
            Where X = (x - Wx) // Sx + 1
        Notes
        -----
        In general, given
          `out` = sliding_window_view(arr,
                                      window_shape=[Wx, (...), Wz],
                                      steps=[Sx, (...), Sz])
           out[ix, (...), iz] = arr[..., ix*Sx:ix*Sx+Wx,  (...), iz*Sz:iz*Sz+Wz]
         Examples
         --------
         >>> import numpy as np
         >>> x = np.arange(9).reshape(3,3)
         >>> x
         array([[0, 1, 2],
                [3, 4, 5],
                [6, 7, 8]])
         >>> y = sliding_window_view(x, window_shape=(2, 2), steps=(1, 1))
         >>> y
         array([[[[0, 1],
                  [3, 4]],
                 [[1, 2],
                  [4, 5]]],
                [[[3, 4],
                  [6, 7]],
                 [[4, 5],
                  [7, 8]]]])
        >>> np.shares_memory(x, y)
         True
        # Performing a neural net style 2D conv (correlation)
        # placing a 4x4 filter with stride-1
        >>> data = np.random.rand(10, 3, 16, 16)  # (N, C, H, W)
        >>> filters = np.random.rand(5, 3, 4, 4)  # (F, C, Hf, Wf)
        >>> windowed_data = sliding_window_view(data,
        ...                                     window_shape=(4, 4),
        ...                                     steps=(1, 1))
        >>> conv_out = np.tensordot(filters,
        ...                         windowed_data,
        ...                         axes=[[1,2,3], [3,4,5]])
        # (F, H', W', N) -> (N, F, H', W')
        >>> conv_out = conv_out.transpose([3,0,1,2])
         """
    import numpy as np
    from numpy.lib.stride_tricks import as_strided
    in_shape = np.array(arr.shape[-len(steps):])  # [x, (...), z]
    window_shape = np.array(window_shape)  # [Wx, (...), Wz]
    steps = np.array(steps)  # [Sx, (...), Sz]
    nbytes = arr.strides[-1]  # size (bytes) of an element in `arr`

    # number of per-byte steps to take to fill window
    window_strides = tuple(np.cumprod(arr.shape[:0:-1])[::-1]) + (1,)
    # number of per-byte steps to take to place window
    step_strides = tuple(window_strides[-len(steps):] * steps)
    # number of bytes to step to populate sliding window view
    strides = tuple(int(i) * nbytes for i in step_strides + window_strides)

    outshape = tuple((in_shape - window_shape) // steps + 1)
    # outshape: ([X, (...), Z], ..., [Wx, (...), Wz])
    outshape = outshape + arr.shape[:-len(steps)] + tuple(window_shape)
    return as_strided(arr, shape=outshape, strides=strides, writeable=False)