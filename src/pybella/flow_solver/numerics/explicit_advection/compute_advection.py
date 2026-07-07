import numba as nb

from ....backends import is_jax_backend
from ....utils import axes
from ....utils import options as opts
from ....utils.slices import get_neighbor_indices
from ...utils.boundary import cell_boundary as bdry_c
from . import advective_flux, recovery, riemann_solver


def _flip_forward(mem):
    """Flip the solution AND the terrain metric so ghost-cell fills and
    flux kernels inside a sweep see consistently oriented arrays."""
    mem.sol.flip_forward()
    if mem.elem.metric is not None:
        mem.elem.metric.flip_forward()


def _flip_backward(mem):
    mem.sol.flip_backward()
    if mem.elem.metric is not None:
        mem.elem.metric.flip_backward()


def strange_splitting(mem, ud, dt, odd, label, writer=None):
    """
    Concise implementation of Strang-splitting advection.
    This function updates the `Sol` solution container with the advected solution in-place.
    """
    time_step = 0.5 * dt
    diagnostics = [writer, mem.node, label] if writer is not None else None

    # Define sweep configurations: (reverse_order, use_diagnostics)
    sweeps = [(not odd, not odd), (odd, odd)]

    for reverse_order, use_diagnostics in sweeps:
        _perform_dimensional_sweep(
            mem,
            ud,
            time_step,
            reverse=reverse_order,
            diagnostics=diagnostics if use_diagnostics else None,
        )

    bdry_c.set_ghost_cells(mem, ud)


def first_order_runge_kutta(mem, ud, dt):
    """
    Function that runs the advection routine with a first-order Runge-Kutta update.
    This function updates the `Sol` solution container with the advected solution in-place.

    Attention
    ---------
    This function is not usually called unless commented out in the :py:meth:`management.data.time_update` routine.
    """
    time_step = dt
    ndim = mem.elem.ndim
    flux = mem.cache.get_flux_containers(mem.elem)

    # Compute fluxes for all dimensions
    for split in range(ndim):
        lmbda = time_step / mem.elem.dxyz[split]
        _flip_forward(mem)
        if mem.elem.iisc[split] > 1:
            flux[split] = _explicit_step_and_flux(mem, ud, lmbda, split, tag="rk")

    # Cache neighbor indices once
    left_idx, right_idx = get_neighbor_indices(mem.elem.ndim)

    # Apply flux updates for all dimensions
    for dim in range(ndim):
        _apply_dimensional_flux_update(mem, ud, dim, time_step, left_idx, right_idx)

    bdry_c.set_ghost_cells(mem, ud)


def _update_solution_variables(
    sol, flux, lmbda, left_idx, right_idx, variables=None, ooJ=None
):
    """
    Helper function to update solution variables with flux differences.

    ooJ: terrain inverse Jacobian (sweep-oriented); the finite-volume cell
    measure is J * dxi, so metric flux differences are divided by J.
    """
    if variables is None:
        variables = ["rho", "rhou", "rhov", "rhow", "rhoX", "rhoY"]

    for var in variables:
        flux_diff = getattr(flux, var)[left_idx] - getattr(flux, var)[right_idx]
        if ooJ is not None:
            flux_diff = ooJ * flux_diff
        current_val = getattr(sol, var)
        setattr(sol, var, current_val + lmbda * flux_diff)


def _explicit_step_and_flux(mem, ud, lmbda, split_step, tag=None):
    """
    For each advection substep, solve the advection problem. For more details, see :ref:`advection_routine`.
    This function updates the solution `Sol` container in-place if a Strang-splitting is used,
    or returns the `flux` data container if a Runge-Kutta method is used.
    """
    bdry_c.set_ghost_cells(mem, ud, step=split_step)

    flux = mem.cache.get_flux_containers(mem.elem)[split_step]

    flux = _compute_flux_and_recovery(mem, flux, ud, lmbda, split_step, tag)

    # pole axis: kill the conservative flux through the zero-area pole faces
    # so no mass/tracer leaks there (Stage F, F2). Over-pole transport is
    # carried by the longitude sweep.
    if ud.bdry_type[split_step] == opts.BdryType.POLE:
        advective_flux.zero_pole_faces(
            flux, advective_flux._ALL_FLUX, int(mem.elem.igs[split_step])
        )

    # Consider caching neighbor indices
    left_idx, right_idx = get_neighbor_indices(mem.elem.ndim)

    if tag != "rk":
        ooJ = mem.elem.metric.ooJ if mem.elem.metric is not None else None
        _update_solution_variables(mem.sol, flux, lmbda, left_idx, right_idx, ooJ=ooJ)

    if tag == "rk":
        return flux


def _compute_flux_and_recovery(mem, flux, ud, lmbda, split_step, tag=None):
    """
    Helper function to compute flux using gradient recovery and HLL solver.

    Returns:
        flux: Computed flux container
    """
    if is_jax_backend(ud):
        from ....backends.jax_ops import advection as jax_advection

        return jax_advection.compute_flux(mem, flux, ud, lmbda, split_step, tag)

    Lefts, Rights = recovery.compute(mem, flux, ud, lmbda, split_step, tag)

    flux = riemann_solver.hll(mem, flux, Lefts, Rights)

    return flux


def _apply_dimensional_flux_update(mem, ud, dim, time_step, left_idx, right_idx):
    """
    Apply flux update for a specific dimension.
    """
    lmbda = time_step / mem.elem.dxyz[dim]
    _flip_forward(mem)
    flux = mem.cache.get_flux_containers(mem.elem)[dim]

    ooJ = mem.elem.metric.ooJ if mem.elem.metric is not None else None
    _update_solution_variables(mem.sol, flux, lmbda, left_idx, right_idx, ooJ=ooJ)

    # Handle special case for the vertical axis
    if dim == axes.vertical_axis(ud):
        updt = lmbda * (flux.rhoX[left_idx] - flux.rhoX[right_idx])
        setattr(mem.sol, "pwchi", updt)


def _perform_dimensional_sweep(mem, ud, time_step, reverse=False, diagnostics=None):
    """
    Perform a dimensional sweep in either forward or reverse order.

    """
    elem, Sol = mem.elem, mem.sol
    ndim = elem.ndim

    # Determine dimension order
    dim_range = range(ndim - 1, -1, -1) if reverse else range(ndim)

    for split in dim_range:
        lmbda = time_step / elem.dxyz[split]

        # Handle solution flipping based on sweep direction
        if reverse:
            if elem.iisc[split] > 1:
                _explicit_step_and_flux(mem, ud, lmbda, split, diagnostics)
            _flip_backward(mem)
        else:
            _flip_forward(mem)
            if elem.iisc[split] > 1:
                _explicit_step_and_flux(mem, ud, lmbda, split, diagnostics)
