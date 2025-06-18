import numpy as np
from ....utils import options as opts
from .common import get_ghost_padding

def set_ghost_nodes(p, node, ud, igs=None):
    if igs is None:
        igs = node.igs
    for dim in range(node.ndim):
        ghost_padding, idx = get_ghost_padding(node.ndim, dim, igs)

        if ud.bdry_type[dim] == opts.BdryType.PERIODIC:
            p[...] = np.pad(p[idx], ghost_padding, periodic_plus_one)
        else:  # ud.bdry_type[dim] == opts.BdryType.WALL:
            p[...] = np.pad(p[idx], ghost_padding, "reflect")

    # if periodic_plus_one
    if node.iicy == 2:  # implying horizontal slices
        pn = p[:, 2, :]
        pn = np.expand_dims(pn, axis=1)
        p[...] = np.repeat(pn, node.icy, axis=1)


def periodic_plus_one(vector, pad_width, iaxis, kwargs=None):
    """
    Taken from the reference:

    Parameters
    ----------
    vector : ndarray
        A rank 1 array already padded with zeros. Padded values are vector `[:iaxis_pad_width[0]] and vector[-iaxis_pad_width[1]:]`.
    iaxis_pad_width : tuple
        A 2-tuple of ints, `iaxis_pad_width[0]` represents the number of values padded at the beginning of vector where `iaxis_pad_width[1]` represents the number of values padded at the end of vector.
    iaxis : int
        The axis currently being calculated.
    kwargs : dict
        Any keyword arguments the function requires.
    References
    ----------
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.pad.html

    """
    if all(pad_width) > 0:
        vector[: pad_width[0] + 1], vector[-pad_width[1] - 1 :] = (
            vector[-pad_width[1] - pad_width[1] - 1 : -pad_width[1]],
            vector[pad_width[0] : pad_width[0] + pad_width[0] + 1].copy(),
        )
    return vector


