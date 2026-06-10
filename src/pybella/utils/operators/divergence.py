import numba as nb
from .. import options as opts
from . import finite_difference


@nb.njit(cache=True)
def compute_2d(u_field, v_field, dx, dy):
    """
    Compute 2D divergence: ∇·F = ∂u/∂x + ∂v/∂y

    Parameters
    ----------
    u_field : np.ndarray
        Field component in x-direction
    v_field : np.ndarray
        Field component in y-direction
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction

    Returns
    -------
    np.ndarray
        Divergence field averaged to cell centers
    """
    # X-direction: ∂u/∂x
    div_x = finite_difference.do_1d(u_field, dx, axis=0)
    # Average to y-cell centers
    div_x = 0.5 * (div_x[:, :-1] + div_x[:, 1:])

    # Y-direction: ∂v/∂y
    div_y = finite_difference.do_1d(v_field, dy, axis=1)
    # Average to x-cell centers
    div_y = 0.5 * (div_y[:-1, :] + div_y[1:, :])

    return div_x + div_y


@nb.njit(cache=True)
def compute_3d_components(u_field, v_field, w_field, dx, dy, dz):
    """
    Compute 3D divergence: ∇·F = ∂u/∂x + ∂v/∂y + ∂w/∂z

    Parameters
    ----------
    u_field : np.ndarray
        Field component in x-direction
    v_field : np.ndarray
        Field component in y-direction
    w_field : np.ndarray
        Field component in z-direction
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
    dz : float
        Grid spacing in z-direction

    Returns
    -------
    tuple
        (div_x, div_y, div_z) - Individual divergence components
    """
    # X-direction: ∂u/∂x
    div_x = finite_difference.do_1d(u_field, dx, axis=0)
    # Average to y-cell centers, then to z-faces
    div_x = 0.5 * (div_x[:, :-1, :] + div_x[:, 1:, :])
    # the legacy "-0.5" here was a sign error (introduced Oct 2021, archive
    # commit 3661b9d); the divergence must be sign-symmetric in all dims
    div_x = 0.5 * (div_x[:, :, :-1] + div_x[:, :, 1:])

    # Y-direction: ∂v/∂y
    div_y = finite_difference.do_1d(v_field, dy, axis=1)
    # Average to x-cell centers, then to z-faces
    div_y = 0.5 * (div_y[:-1, :, :] + div_y[1:, :, :])
    div_y = 0.5 * (div_y[:, :, :-1] + div_y[:, :, 1:])

    # Z-direction: ∂w/∂z
    div_z = finite_difference.do_1d(w_field, dz, axis=2)
    # Average to cell centers
    div_z = 0.5 * (div_z[:-1, :, :] + div_z[1:, :, :])
    div_z = 0.5 * (div_z[:, :-1, :] + div_z[:, 1:, :])

    return div_x, div_y, div_z


@nb.njit(cache=True)
def compute_3d_sum(u_field, v_field, w_field, dx, dy, dz):
    """
    Compute total 3D divergence: ∇·F = ∂u/∂x + ∂v/∂y + ∂w/∂z

    Parameters
    ----------
    u_field : np.ndarray
        Field component in x-direction
    v_field : np.ndarray
        Field component in y-direction
    w_field : np.ndarray
        Field component in z-direction
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
    dz : float
        Grid spacing in z-direction

    Returns
    -------
    np.ndarray
        Total divergence field
    """
    div_x, div_y, div_z = compute_3d_components(u_field, v_field, w_field, dx, dy, dz)
    return div_x + div_y + div_z


def compute_at_nodes(rhs, elem, sol, ud):
    """Main divergence function - handles boundary conditions and calls JIT-compiled core."""
    ndim = elem.ndim

    # Handle boundary conditions
    if not hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        if (
            ud.bdry_type[1] == opts.BdryType.WALL
            or ud.bdry_type[1] == opts.BdryType.RAYLEIGH
        ):
            sol.rhou[:, :2, ...] = 0.0
            sol.rhov[:, :2, ...] = 0.0
            sol.rhow[:, :2, ...] = 0.0
            sol.rhou[:, -2:, ...] = 0.0
            sol.rhov[:, -2:, ...] = 0.0
            sol.rhow[:, -2:, ...] = 0.0

    # Call appropriate JIT-compiled function
    if ndim == 2:
        rhs[:] = _momentum_pot_temp_divergence_2d_jit(
            sol.rho, sol.rhou, sol.rhov, sol.rhoY, elem.dx, elem.dy
        )
    else:
        _momentum_pot_temp_divergence_3d_jit(
            rhs,
            sol.rho,
            sol.rhou,
            sol.rhov,
            sol.rhow,
            sol.rhoY,
            elem.dx,
            elem.dy,
            elem.dz,
        )

    return rhs


@nb.njit(cache=True)
def _momentum_pot_temp_divergence_2d_jit(rho, rhou, rhov, rhoY, dx, dy):
    """
    JIT-compiled 2D momentum-potential temperature divergence calculation.
    Computes ∇·(ρu θ, ρv θ) where θ = ρY/ρ is the potential temperature.
    """
    # Calculate potential temperature θ = ρY / ρ
    theta = rhoY / rho

    # Compute momentum-potential temperature flux components
    rhou_theta = rhou * theta  # x-momentum flux weighted by potential temperature
    rhov_theta = rhov * theta  # y-momentum flux weighted by potential temperature

    # Use generic divergence operator
    return compute_2d(rhou_theta, rhov_theta, dx, dy)


@nb.njit(cache=True)
def _momentum_pot_temp_divergence_3d_jit(rhs, rho, rhou, rhov, rhow, rhoY, dx, dy, dz):
    """
    JIT-compiled 3D momentum-potential temperature divergence calculation.
    Computes ∇·(ρu θ, ρv θ, ρw θ) where θ = ρY/ρ is the potential temperature.
    """
    # Calculate potential temperature θ = ρY / ρ
    theta = rhoY / rho

    # Compute momentum-potential temperature flux components
    rhou_theta = rhou * theta  # x-momentum flux weighted by potential temperature
    rhov_theta = rhov * theta  # y-momentum flux weighted by potential temperature
    rhow_theta = rhow * theta  # z-momentum flux weighted by potential temperature

    # Use generic total divergence operator; rhs is interior-sized (node.isc),
    # which is exactly the shape the cell-array differences produce
    rhs[:, :, :] = compute_3d_sum(rhou_theta, rhov_theta, rhow_theta, dx, dy, dz)
