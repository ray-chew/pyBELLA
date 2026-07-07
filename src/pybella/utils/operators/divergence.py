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
            if ud.bdry_type[dim] in (
                opts.BdryType.WALL,
                opts.BdryType.RAYLEIGH,
                opts.BdryType.POLE,
            ):
                # POLE is one-sided like a wall for the elliptic rhs (Stage
                # F): zero the beyond-pole ghost momenta so the pole-node
                # divergence matches the one-sided operator; the pole-ring
                # collapse then sums these into the master equation
                lo, hi = axes.wall_slabs(ndim, dim)
                for field in (sol.rhou, sol.rhov, sol.rhow):
                    field[lo] = 0.0
                    field[hi] = 0.0

    # Call appropriate JIT-compiled function
    if ndim == 2:
        if elem.metric is not None:
            # terrain: general curvilinear fluxes F_a = N_a . (theta m);
            # for the vertical-line metric this is bit-exactly the legacy
            # J-weighted / contravariant construction (Phase-0 contract,
            # test_scripts/test_metric_reduction.py)
            m = elem.metric
            f_h1, f_v = _normal_fluxes_2d_jit(
                sol.rho,
                sol.rhoY,
                sol.rhou,
                sol.rhov,
                m.N[0][m.cart_v],
                m.N[0][m.cart_haxes[0]],
                m.N[1][m.cart_v],
                m.N[1][m.cart_haxes[0]],
            )
            rhs[:] = compute_2d(f_h1, f_v, elem.dx, elem.dy)
        else:
            rhs[:] = _momentum_pot_temp_divergence_2d_jit(
                sol.rho, sol.rhou, sol.rhov, sol.rhoY, elem.dx, elem.dy
            )
    elif elem.metric is not None:
        # terrain: rhs = J grad.F with the general curvilinear fluxes
        # F_a = N_a . (theta m) per array axis a (Cartesian components,
        # vertical-first contraction); the differencing stencils are
        # unchanged
        m = elem.metric
        moms = (sol.rhou, sol.rhov, sol.rhow)
        cv = m.cart_v
        ch1, ch2 = m.cart_haxes
        flux = [
            _normal_flux_3d_jit(
                sol.rho,
                sol.rhoY,
                moms[cv],
                moms[ch1],
                moms[ch2],
                m.N[a][cv],
                m.N[a][ch1],
                m.N[a][ch2],
            )
            for a in range(3)
        ]
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
def _normal_flux_3d_jit(rho, rhoY, mom_v, mom_h1, mom_h2, n_v, n_h1, n_h2):
    """One general curvilinear flux component F_a = N_a . (theta m).

    ``n_*`` are the Cartesian components of the area normal N_a of the
    xi_a = const surface, passed vertical-first — the contraction order
    that reduces BIT-EXACTLY to the legacy J-weighted / contravariant
    fluxes for vertical-line metrics (Phase-0 contract,
    ``test_scripts/test_metric_reduction.py``). The plain divergence of
    the three F_a equals J grad.F in physical space.
    """
    theta = rhoY / rho
    return n_v * (mom_v * theta) + n_h1 * (mom_h1 * theta) + n_h2 * (mom_h2 * theta)


@nb.njit(cache=True)
def _normal_fluxes_2d_jit(rho, rhoY, mom_h1, mom_v, n1_v, n1_h1, n2_v, n2_h1):
    """2D restriction of :func:`_normal_flux_3d_jit` (both components).

    Vertical-first contraction; reduces bit-exactly to the legacy
    (J f_h1, f_v - G1 f_h1) pair for vertical-line metrics.
    """
    theta = rhoY / rho
    f_h1 = mom_h1 * theta
    f_v = mom_v * theta
    return n1_v * f_v + n1_h1 * f_h1, n2_v * f_v + n2_h1 * f_h1


@nb.njit(cache=True)
def _metric_contravariant_fluxes_jit(rho, rhoY, mom_h1, mom_v, mom_h2, J, G1, G2):
    """Legacy vertical-line flux components (reduction-contract reference).

    Kept as the reference the general :func:`_normal_flux_3d_jit` path is
    pinned against in ``test_scripts/test_metric_reduction.py``; the
    production divergence no longer calls it. Role-ordered inputs/outputs
    (h1, v, h2): J-weighted horizontal fluxes and the contravariant
    vertical flux f_v = theta * (mom_v - G1 mom_h1 - G2 mom_h2).
    """
    theta = rhoY / rho
    f_h1 = mom_h1 * theta
    f_h2 = mom_h2 * theta
    f_v = mom_v * theta - G1 * f_h1 - G2 * f_h2
    return J * f_h1, f_v, J * f_h2


@nb.njit(cache=True)
def _metric_contravariant_fluxes_2d_jit(rho, rhoY, mom_h1, mom_v, J, G1):
    """2D restriction of :func:`_metric_contravariant_fluxes_jit`.

    Reduction-contract reference only, like the 3D variant.
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
