import numpy as np

from ....utils import options as opts
from ....backends import is_jax_backend


def get_ghost_padding(ndim, dim, igs):
    """
    For a given direction, return the number of ghost cells to pad the current direction, and the index slice of the inner array without the ghost cells.

    Parameters
    ----------
    ndim : int
        Number of dimensions for the problem.
    dim : int
        Current dimension to update
    igs : list
        A list of number of ghost cells in all dimensions, e.g. `[2,2,2]` for 2 ghost cells in the x, y, and z directions.

    Returns
    -------
    tuple
        Number of ghost cells in the current dimension at both edges.
    tuple
        Index slice of the inner domain of the array.

    """
    ghost_padding = [(0, 0)] * ndim
    ghost_padding[dim] = (igs[dim], igs[dim])

    padded_idx = np.empty((ndim), dtype=object)
    for idim in range(ndim):
        padded_idx[idim] = slice(igs[idim], -igs[idim])
    padded_idx[dim] = slice(None)

    inner_domain = [slice(None)] * ndim
    inner_domain[dim] = slice(igs[dim], -igs[dim])

    return tuple(ghost_padding), tuple(inner_domain)


def pole_source_indices(ncx, ncz, ig, nodal):
    """Index maps realizing the lat-lon pole fold.

    A ghost cell/node past |phi| = pi/2 covers the physical point on the
    FAR side of the pole, at longitude lambda + pi and the mirror latitude.
    With GLOBAL Cartesian momenta the exchange is a pure index remap (no
    vector rotation): the source is an INTERIOR cell obtained by

    * shifting longitude by half the ring (``+pi`` = N/2 cells, always
      landing inside the interior lambda range — so the fill never reads a
      stale lambda-ghost column, order-independent within a sweep), and
    * mirroring latitude across the pole (cells mirror about the boundary
      FACE, nodes reflect about the pole NODE).

    ``ncx``/``ncz`` are the full padded lambda/phi extents, ``ig`` the ghost
    width. Returns ``(src_lam, src_phi, ghost_phi_slabs)`` where ``src_lam``
    is length ``ncx``, ``src_phi`` length ``ncz`` (identity on interior phi
    rows — only the ghost slabs are overwritten by the caller).
    """
    # interior lambda period in CELLS (nodes carry the +pi seam duplicate,
    # so their interior count is one larger); the +pi remap is N/2 cells
    N = (ncx - 2 * ig) - (1 if nodal else 0)
    j = np.arange(ncx)
    src_lam = ig + ((j - ig + N // 2) % N)

    src_phi = np.arange(ncz)
    lo = np.arange(0, ig)
    hi = np.arange(ncz - ig, ncz)
    if nodal:  # reflect about the pole node
        src_phi[lo] = 2 * ig - lo
        src_phi[hi] = 2 * (ncz - 1 - ig) - hi
    else:  # mirror about the wall face
        src_phi[lo] = 2 * ig - 1 - lo
        src_phi[hi] = 2 * (ncz - ig) - 1 - hi
    return src_lam, src_phi, (slice(0, ig), slice(ncz - ig, ncz))


def pole_exchange_field(arr, lam_axis, phi_axis, src_lam, src_phi, ghost_slabs):
    """In-place pole fold of one field's phi-ghost slabs (see
    :func:`pole_source_indices`). The source is gathered from interior
    lambda/phi only, so no vector components rotate and no sign flips."""
    gathered = np.take(np.take(arr, src_lam, axis=lam_axis), src_phi, axis=phi_axis)
    for slab in ghost_slabs:
        dst = [slice(None)] * arr.ndim
        dst[phi_axis] = slab
        arr[tuple(dst)] = gathered[tuple(dst)]


def scale_wall_node_values(rhs, node, ud, factor=0.5):
    """Scale values at wall boundary nodes by a given factor."""
    if is_jax_backend(ud):
        from ....backends.jax_ops import boundary as jax_boundary

        return jax_boundary.scale_wall_node_values(rhs, node, ud, factor=factor)

    ndim = node.ndim
    igs = node.igs

    for dim in range(ndim):
        # Check if this dimension has wall boundaries
        is_wall = (
            ud.bdry_type[dim] == opts.BdryType.WALL
            or ud.bdry_type[dim] == opts.BdryType.RAYLEIGH
        )

        if is_wall:
            # Create index for all dimensions
            idx = [slice(igs[d], -igs[d]) for d in range(ndim)]

            # Scale first and last interior nodes in this dimension
            for boundary_idx in [igs[dim], -igs[dim] - 1]:
                idx[dim] = boundary_idx
                rhs[tuple(idx)] *= factor
                # rhs = rhs.at[tuple(idx)].multiply(factor)

    return rhs
