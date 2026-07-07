import numpy as np
from ....utils import axes
from ....utils import options as opts
from . import common as bdry
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
        elif ud.bdry_type[dim] == opts.BdryType.POLE:
            _apply_pole_nodes(p, node, dim)
        else:  # ud.bdry_type[dim] == opts.BdryType.WALL:
            p[...] = np.pad(p[idx], ghost_padding, "reflect")

    # quasi-2D: broadcast the single interior layer across any degenerate
    # axis (historically hardcoded to axis 1 / iicy == 2)
    for dim in axes.degenerate_axes(node):
        slc = [slice(None)] * p.ndim
        slc[dim] = node.igs[dim]
        pn = np.expand_dims(p[tuple(slc)], axis=dim)
        p[...] = np.repeat(pn, node.sc[dim], axis=dim)


def _apply_pole_nodes(p, node, dim):
    """Pole node exchange + ring-consistency on the phi axis (Stage F, F1).

    Same pure index remap as the cell fill, on the node pressure field.
    Additionally forces the two pole-node rows (phi = +-pi/2, one physical
    point per radius) to their lambda-ring MEAN, keeping p2 single-valued
    at the pole between elliptic solves (the collapse in F4 produces it
    single-valued; this holds it so between solves). Nodes are never
    sweep-flipped, so lambda is array axis 0 and phi is ``dim``.
    """
    lam_axis, phi_axis = 0, dim
    ig = int(node.igs[0])
    ncx = p.shape[lam_axis]
    ncz = p.shape[phi_axis]
    src_lam, src_phi, slabs = bdry.pole_source_indices(ncx, ncz, ig, nodal=True)
    bdry.pole_exchange_field(p, lam_axis, phi_axis, src_lam, src_phi, slabs)

    # ring-consistency: each pole row -> its unique-longitude ring mean
    # (exclude the +pi seam duplicate node); broadcast over all longitudes
    n_int = (ncx - 2 * ig) - 1  # interior lambda intervals = unique nodes
    for pole_row in (ig, ncz - 1 - ig):
        src = [slice(None)] * p.ndim
        src[phi_axis] = slice(pole_row, pole_row + 1)
        src[lam_axis] = slice(ig, ig + n_int)
        mean = p[tuple(src)].mean(axis=lam_axis, keepdims=True)
        dst = [slice(None)] * p.ndim
        dst[phi_axis] = slice(pole_row, pole_row + 1)
        p[tuple(dst)] = mean


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
