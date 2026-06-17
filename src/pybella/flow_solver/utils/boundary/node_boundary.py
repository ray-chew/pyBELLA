import numpy as np
from ....utils import axes
from ....utils import options as opts
from .common import get_ghost_padding
from ....backends import is_jax_backend


def set_ghost_nodes(p, node, ud, igs=None):
    if is_jax_backend(ud):
        from ....backends.jax_ops import boundary as jax_boundary

        return jax_boundary.set_ghost_nodes(p, node, ud, igs=igs)

    if igs is None:
        igs = node.igs
    for dim in range(node.ndim):
        ghost_padding, idx = get_ghost_padding(node.ndim, dim, igs)

        if ud.bdry_type[dim] == opts.BdryType.PERIODIC:
            p[...] = np.pad(p[idx], ghost_padding, periodic_plus_one)
        else:  # ud.bdry_type[dim] == opts.BdryType.WALL:
            p[...] = np.pad(p[idx], ghost_padding, "reflect")

    # quasi-2D: broadcast the single interior layer across any degenerate
    # axis (historically hardcoded to axis 1 / iicy == 2)
    for dim in axes.degenerate_axes(node):
        slc = [slice(None)] * p.ndim
        slc[dim] = node.igs[dim]
        pn = np.expand_dims(p[tuple(slc)], axis=dim)
        p[...] = np.repeat(pn, node.sc[dim], axis=dim)


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
