from ....utils.slices import get_neighbor_indices
from ...utils import boundary as bdry
from . import recovery, riemann_solver


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

    bdry.set_explicit_boundary_data(mem.sol, mem.elem, ud, mem.th, mem.mpv)


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

    # Compute fluxes for all dimensions
    for split in range(ndim):
        lmbda = time_step / mem.elem.dxyz[split]
        mem.sol.flip_forward()
        if mem.elem.iisc[split] > 1:
            mem.flux[split] = _explicit_step_and_flux(mem, ud, lmbda, split, tag="rk")

    # Cache neighbor indices once
    left_idx, right_idx = get_neighbor_indices(mem.elem.ndim)

    # Apply flux updates for all dimensions
    for dim in range(ndim):
        _apply_dimensional_flux_update(mem, dim, time_step, left_idx, right_idx)

    bdry.set_explicit_boundary_data(mem.sol, mem.elem, ud, mem.th, mem.mpv)


def _update_solution_variables(sol, flux, lmbda, left_idx, right_idx, variables=None):
    """
    Helper function to update solution variables with flux differences.

    """
    if variables is None:
        variables = ["rho", "rhou", "rhov", "rhow", "rhoX", "rhoY"]

    for var in variables:
        flux_diff = getattr(flux, var)[left_idx] - getattr(flux, var)[right_idx]
        current_val = getattr(sol, var)
        setattr(sol, var, current_val + lmbda * flux_diff)


def _explicit_step_and_flux(mem, ud, lmbda, split_step, tag=None):
    """
    For each advection substep, solve the advection problem. For more details, see :ref:`advection_routine`.
    This function updates the solution `Sol` container in-place if a Strang-splitting is used,
    or returns the `flux` data container if a Runge-Kutta method is used.
    """
    flux = _compute_flux_and_recovery(mem, ud, lmbda, split_step, tag)

    # Cache neighbor indices (consider moving this to initialization if called frequently)
    left_idx, right_idx = get_neighbor_indices(mem.elem.ndim)

    if tag != "rk":
        _update_solution_variables(mem.sol, flux, lmbda, left_idx, right_idx)

    bdry.set_explicit_boundary_data(
        mem.sol, mem.elem, ud, mem.th, mem.mpv, step=split_step
    )

    if tag == "rk":
        return flux


def _compute_flux_and_recovery(mem, ud, lmbda, split_step, tag=None):
    """
    Helper function to compute flux using gradient recovery and HLL solver.

    Returns:
        flux: Computed flux container
    """
    flux = mem.flux[split_step]

    bdry.set_explicit_boundary_data(
        mem.sol, mem.elem, ud, mem.th, mem.mpv, step=split_step
    )

    Lefts, Rights = recovery.compute(mem, ud, lmbda, split_step, tag)

    flux = riemann_solver.hll(mem, flux, Lefts, Rights)

    return flux


def _apply_dimensional_flux_update(mem, dim, time_step, left_idx, right_idx):
    """
    Apply flux update for a specific dimension.
    """
    lmbda = time_step / mem.elem.dxyz[dim]
    mem.sol.flip_forward()

    _update_solution_variables(mem.sol, mem.flux[dim], lmbda, left_idx, right_idx)

    # Handle special case for vertical axis
    if dim == 1:
        updt = lmbda * (mem.flux[dim].rhoX[left_idx] - mem.flux[dim].rhoX[right_idx])
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
            Sol.flip_backward()
        else:
            Sol.flip_forward()
            if elem.iisc[split] > 1:
                _explicit_step_and_flux(mem, ud, lmbda, split, diagnostics)
