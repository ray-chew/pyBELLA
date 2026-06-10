import numba as nb
from .. import axes
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

    # Handle boundary conditions: zero the momenta in the two boundary
    # slabs of every WALL/RAYLEIGH axis (historically vertical-only, which
    # left the x-WALL elliptic path broken). The slabs are the ghost
    # layers, so with terrain this also zeroes the ghost contravariant
    # fluxes (they are formed from the momenta) — the same wall treatment
    # as the uniform-Cartesian path.
    if not hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        for dim in range(ndim):
            if (
                ud.bdry_type[dim] == opts.BdryType.WALL
                or ud.bdry_type[dim] == opts.BdryType.RAYLEIGH
            ):
                lo, hi = axes.wall_slabs(ndim, dim)
                for field in (sol.rhou, sol.rhov, sol.rhow):
                    field[lo] = 0.0
                    field[hi] = 0.0

    # Call appropriate JIT-compiled function
    if ndim == 2:
        if elem.metric is not None:
            # terrain: same contravariant/J-weighted construction as the 3D
            # branch below, minus the second-horizontal leg (haxes = (0, None),
            # vaxis = 1 — enforced by axes.validate in 2D)
            m = elem.metric
            f_h1, f_v = _metric_contravariant_fluxes_2d_jit(
                sol.rho, sol.rhoY, sol.rhou, sol.rhov, m.J, m.G1
            )
            rhs[:] = compute_2d(f_h1, f_v, elem.dx, elem.dy)
        else:
            rhs[:] = _momentum_pot_temp_divergence_2d_jit(
                sol.rho, sol.rhou, sol.rhov, sol.rhoY, elem.dx, elem.dy
            )
    elif elem.metric is not None:
        # terrain: rhs = J grad.F with J-weighted horizontal fluxes and the
        # contravariant vertical flux F_v - G1 F_h1 - G2 F_h2 (role space);
        # the differencing stencils are unchanged
        m = elem.metric
        moms = (sol.rhou, sol.rhov, sol.rhow)
        a_h1, a_h2 = m.haxes
        f_h1, f_v, f_h2 = _metric_contravariant_fluxes_jit(
            sol.rho,
            sol.rhoY,
            moms[a_h1],
            moms[m.vaxis],
            moms[a_h2],
            m.J,
            m.G1,
            m.G2,
        )
        flux = [None, None, None]
        flux[a_h1], flux[m.vaxis], flux[a_h2] = f_h1, f_v, f_h2
        rhs[:, :, :] = compute_3d_sum(
            flux[0], flux[1], flux[2], elem.dx, elem.dy, elem.dz
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
def _metric_contravariant_fluxes_jit(rho, rhoY, mom_h1, mom_v, mom_h2, J, G1, G2):
    """Terrain-following flux components of the theta-weighted momentum.

    Role-ordered inputs/outputs (h1, v, h2). Returns the J-weighted
    horizontal fluxes and the contravariant vertical flux

        f_v = theta * (mom_v - G1 mom_h1 - G2 mom_h2)

    such that the plain divergence of (J f_h1, f_v, J f_h2) equals
    J grad.F in physical space.
    """
    theta = rhoY / rho
    f_h1 = mom_h1 * theta
    f_h2 = mom_h2 * theta
    f_v = mom_v * theta - G1 * f_h1 - G2 * f_h2
    return J * f_h1, f_v, J * f_h2


@nb.njit(cache=True)
def _metric_contravariant_fluxes_2d_jit(rho, rhoY, mom_h1, mom_v, J, G1):
    """2D restriction of :func:`_metric_contravariant_fluxes_jit`.

    Returns the J-weighted horizontal flux and the contravariant vertical
    flux f_v = theta * (mom_v - G1 mom_h1) such that the plain 2D
    divergence of (J f_h1, f_v) equals J grad.F in physical space.
    """
    theta = rhoY / rho
    f_h1 = mom_h1 * theta
    f_v = mom_v * theta - G1 * f_h1
    return J * f_h1, f_v


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
