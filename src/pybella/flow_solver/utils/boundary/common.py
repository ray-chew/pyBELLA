import numpy as np

from ....utils import options as opts


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


def scale_wall_node_values(rhs, node, ud, factor=0.5):
    """Scale values at wall boundary nodes by a given factor."""
    if getattr(ud, "backend", "numpy") in ("jax", "jax-device"):
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
