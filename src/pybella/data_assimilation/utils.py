import copy

import numpy as np

from ..utils import options as opts
from . import ensemble_access


def ensemble_inflation(results, attributes, factor, N):
    for attribute in attributes:
        mean = ensemble_access.stack_fields(results, attribute)
        mean = np.mean(mean, axis=0)
        for n in range(N):
            inflation = mean + factor * (
                ensemble_access.get_field(results[n], attribute) - mean
            )
            ensemble_access.set_field(results[n], attribute, inflation)


def boundary_mask(ud, elem, node, pad_X, pad_Y):
    """
    Returns a mask for the underlying cellular and nodal grids padded such that the size of the local subdomain has been accounted for. For ghost cells on wall boundaries, values are 0.0 and 1.0 for periodic.
    """

    pads = [pad_X, pad_Y]
    cmask = np.ones(elem.iisc).squeeze()
    nmask = np.ones(node.iisc).squeeze()

    for dim in range(elem.ndim):
        ghost_padding = [[0, 0]] * elem.ndim
        ghost_padding[dim] = [pads[dim], pads[dim]]

        if ud.bdry_type[dim] == opts.BdryType.PERIODIC:
            cmask = np.pad(cmask, ghost_padding, mode="constant", constant_values=(1.0))
            nmask = np.pad(nmask, ghost_padding, mode="constant", constant_values=(1.0))

        elif ud.bdry_type[dim] == opts.BdryType.WALL:
            cmask = np.pad(cmask, ghost_padding, mode="constant", constant_values=(0.0))
            nmask = np.pad(nmask, ghost_padding, mode="constant", constant_values=(0.0))

    return cmask.astype("bool"), nmask.astype("bool")


# ref: https://gist.github.com/meowklaski/4bda7c86c6168f3557657d5fb0b5395a
def sliding_window_view(arr, window_shape, steps):
    """
    Produce a view from a sliding, striding window over `arr`.
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
    view of `arr`, shape=([X, (...), Z], ..., [Wx, (...), Wz]), where X = (x - Wx) // Sx + 1

    Note
    -----
    In general, given::

        out = sliding_window_view(arr,
                                    window_shape=[Wx, (...), Wz],
                                    steps=[Sx, (...), Sz])
        out[ix, (...), iz] = arr[..., ix*Sx:ix*Sx+Wx,  (...), iz*Sz:iz*Sz+Wz]

    This function is taken from:
    https://gist.github.com/meowklaski/4bda7c86c6168f3557657d5fb0b5395a

    Example
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

    from numpy.lib.stride_tricks import as_strided

    in_shape = np.array(arr.shape[-len(steps) :])  # [x, (...), z]
    window_shape = np.array(window_shape)  # [Wx, (...), Wz]
    steps = np.array(steps)  # [Sx, (...), Sz]
    nbytes = arr.strides[-1]  # size (bytes) of an element in `arr`

    # number of per-byte steps to take to fill window
    window_strides = tuple(np.cumprod(arr.shape[:0:-1])[::-1]) + (1,)
    # number of per-byte steps to take to place window
    step_strides = tuple(window_strides[-len(steps) :] * steps)
    # number of bytes to step to populate sliding window view
    strides = tuple(int(i) * nbytes for i in step_strides + window_strides)

    outshape = tuple((in_shape - window_shape) // steps + 1)
    # outshape: ([X, (...), Z], ..., [Wx, (...), Wz])
    outshape = outshape + arr.shape[: -len(steps)] + tuple(window_shape)
    return as_strided(arr, shape=outshape, strides=strides, writeable=False)


def sparse_obs_selector(obs, elem, node, ud, dap):
    sparse_obs = dap.sparse_obs

    if not sparse_obs or len(dap.da_times) == 0:
        mask = copy.deepcopy(obs)
        for tt, mask_t in enumerate(mask):
            for key, _ in mask_t.items():
                mask[tt][key][...] = 0.0
        return obs, mask

    else:
        sparse_obs_by_attr = dap.sparse_obs_by_attr
        seeds = dap.sparse_obs_seeds
        K = dap.obs_frac

        # define inner and outer domains in 2D
        i2 = (slice(elem.igx, -elem.igx), slice(elem.igy, -elem.igy))
        i0 = (
            slice(
                None,
            ),
            slice(
                None,
            ),
        )

        # get inner domain size
        Ncx, Ncy = elem.iicx, elem.iicy
        Nnx, Nny = node.iicx, node.iicy

        Nc = Ncx * Ncy
        Nn = Nnx * Nny
        Kc = Nc * K
        Kn = Nn * K
        Kc = int(np.ceil(Kc))
        Kn = int(np.ceil(Kn))

        Xc, Yc = elem.x, elem.y
        Xc, Yc = np.meshgrid(Xc, Yc)
        Xn, Yn = node.x, node.y
        Xn, Yn = np.meshgrid(Xn, Yn)

        mask_arr = copy.deepcopy(obs)
        # obs_noisy_interp = deepcopy(obs)

        # obs is a list of dictionaries, list length da_len, dictionary length attr_len.
        for tt, obs_t in enumerate(obs):
            attr_cnt = tt
            for key, value in obs_t.items():
                if key == "p2_nodes":
                    grid_x, grid_y = Xn[i0], Yn[i0]
                    K, N = Kn, Nn
                    Nx, Ny = Nnx, Nny
                else:
                    grid_x, grid_y = Xc[i0], Yc[i0]
                    K, N = Kc, Nc
                    Nx, Ny = Ncx, Ncy
                grid_x, grid_y = grid_x.T, grid_y.T
                np.random.seed(seeds[attr_cnt])

                # ref: method for generating random boolean array given a probability of 1's and 0's.
                # https://stackoverflow.com/questions/19597473/binary-random-array-with-a-specific-proportion-of-ones/19597805

                # append mask array with new seed
                mask = np.array([0] * K + [1] * (N - K))
                np.random.shuffle(mask)
                mask = mask.reshape(Nx, Ny)
                mask = np.pad(
                    mask,
                    ((elem.igx, elem.igx), (elem.igy, elem.igy)),
                    mode="constant",
                    constant_values=0.0,
                )

                # values = np.ma.array(value[i0], mask=mask).compressed()
                # X = np.ma.array(grid_x, mask=mask).compressed()
                # Y = np.ma.array(grid_y, mask=mask).compressed()

                # points = np.zeros((len(values),2))
                # points[:,0] = X[...].flatten()
                # points[:,1] = Y[...].flatten()

                # values = griddata(points, values, (grid_x, grid_y), method='cubic')

                # if dap.obs_frac < 1.0:
                # obs_noisy_interp[tt][key][...] = 0.0
                # obs_noisy_interp[tt][key][i0] = values

                mask_arr[tt][key][...] = 1
                mask_arr[tt][key][i2] = mask[i2]

                if sparse_obs_by_attr:
                    attr_cnt += 1

        for mask_at_t in mask_arr:
            for key, mask in mask_at_t.items():
                if key == "p2_nodes":
                    Nx, Ny = Nnx, Nny
                else:
                    Nx, Ny = Ncx, Ncy
                assert (mask.shape[0] * mask.shape[1]) - mask.sum() == np.ceil(
                    Nx * Ny * dap.obs_frac
                ), "Mask sparsity does not match obs_frac defined"

        # return obs_noisy_interp, mask_arr
        return mask_arr


def obs_noiser(obs, mask, dap, rloc, elem):
    if dap.add_obs_noise:
        assert isinstance(dap.obs_noise, dict), "obs_noise has to be dict"
        assert len(dap.obs_noise) == len(
            dap.obs_attributes
        ), "obs_noise length has to be equal to obs_attributes len"

        obs_covar_c = np.zeros((len(dap.da_times), rloc.cattr_len))
        obs_covar_n = np.zeros((len(dap.da_times), rloc.nattr_len))
        obs_noisy = copy.deepcopy(obs)

        std_dev = np.zeros((obs.shape[0], len(dap.obs_noise_seeds)))

        # define inner domain in 2D VS
        i2 = (slice(elem.igx, -elem.igx), slice(elem.igy, -elem.igy))

        for tt, obs_t in enumerate(obs):
            attr_cnt = 0
            for key, value in obs_t.items():
                # Get inner 2D domain of sparse obs field
                value = value[i2]
                mask_at_t = mask[tt][key][i2]
                value = np.ma.array(value, mask=mask_at_t)

                seed = dap.obs_noise_seeds[attr_cnt]
                np.random.seed(seed)

                if dap.noise_type == "AmpCov":
                    field_var = dap.obs_noise[key] * np.abs(value.max() - value.min())
                    field_sd = field_var**0.5
                elif dap.noise_type == "VarCov":
                    field_var = (
                        dap.obs_noise[key] * ((value - value.mean()) ** 2).mean()
                    )
                    field_sd = field_var**0.5

                elif dap.noise_type == "FixCov":
                    field_sd = 1.0 * dap.obs_noise[key]

                # here, we take the fraction defined by obs_noise multiplied by the maximum value of the observation as the standard deviation of the measurement noise.

                std_dev[tt, attr_cnt] = field_sd

                attr_cnt += 1

                # var = std_dev**2

        if dap.noise_type == "VarCov" or dap.noise_type == "AmpCov":
            sd = std_dev.mean(axis=0, keepdims=True)
            std_dev[:, ...] = sd

        for tt, obs_t in enumerate(obs):
            ccnt, ncnt, attr_cnt = 0, 0, 0
            for key, value in obs_t.items():
                shp = value.shape
                sd = std_dev[tt, attr_cnt]

                # print(tt, key, sd**2, np.abs(value.max()-value.min()))

                # generate gaussian noise for observations.
                noise = np.random.normal(0.0, sd, size=(shp))

                # add noise onto observation
                obs_noisy[tt][key][...] += noise

                var = sd**2

                if var == 0:
                    var = -np.inf
                if key in rloc.ca:
                    obs_covar_c[tt, ccnt] = var
                    ccnt += 1
                else:
                    obs_covar_n[tt, ncnt] = var
                    ncnt += 1
                attr_cnt += 1

        obs_covar = [obs_covar_c, obs_covar_n]
        return obs_noisy, obs_covar

    else:
        return obs, None
