"""JAX twin of :mod:`pybella.utils.operators.divergence`.

The numpy ``compute_at_nodes(rhs, elem, sol, ud)`` mutates ``rhs`` and the
wall slabs of the momenta in place. The JAX twin is functional: it takes
plain arrays plus static flags and returns ``(rhs, momenta)`` with the
slab-zeroed momenta, leaving the inputs untouched. Use
:func:`wall_zero_dims` to derive the static flags from ``ud`` host-side.
"""

import jax.numpy as jnp

from pybella.utils import axes
from pybella.utils import options as opts

from . import finite_difference


def compute_2d(u_field, v_field, dx, dy):
    """Compute 2D divergence averaged to cell centers."""
    div_x = finite_difference.do_1d(u_field, dx, axis=0)
    div_x = 0.5 * (div_x[:, :-1] + div_x[:, 1:])

    div_y = finite_difference.do_1d(v_field, dy, axis=1)
    div_y = 0.5 * (div_y[:-1, :] + div_y[1:, :])

    return div_x + div_y


def compute_3d_components(u_field, v_field, w_field, dx, dy, dz):
    """Compute the three 3D divergence components, averaged transversally."""
    div_x = finite_difference.do_1d(u_field, dx, axis=0)
    div_x = 0.5 * (div_x[:, :-1, :] + div_x[:, 1:, :])
    div_x = 0.5 * (div_x[:, :, :-1] + div_x[:, :, 1:])

    div_y = finite_difference.do_1d(v_field, dy, axis=1)
    div_y = 0.5 * (div_y[:-1, :, :] + div_y[1:, :, :])
    div_y = 0.5 * (div_y[:, :, :-1] + div_y[:, :, 1:])

    div_z = finite_difference.do_1d(w_field, dz, axis=2)
    div_z = 0.5 * (div_z[:-1, :, :] + div_z[1:, :, :])
    div_z = 0.5 * (div_z[:, :-1, :] + div_z[:, 1:, :])

    return div_x, div_y, div_z


def compute_3d_sum(u_field, v_field, w_field, dx, dy, dz):
    """Compute total 3D divergence."""
    div_x, div_y, div_z = compute_3d_components(u_field, v_field, w_field, dx, dy, dz)
    return div_x + div_y + div_z


def _normal_flux_3d(rho, rhoY, mom_v, mom_h1, mom_h2, n_v, n_h1, n_h2):
    """General curvilinear flux F_a = N_a . (theta m), vertical-first.

    Twin of the numpy ``_normal_flux_3d_jit``: same contraction order so
    the vertical-line reduction stays bit-exact on both backends.
    """
    theta = rhoY / rho
    return n_v * (mom_v * theta) + n_h1 * (mom_h1 * theta) + n_h2 * (mom_h2 * theta)


def _normal_fluxes_2d(rho, rhoY, mom_h1, mom_v, n1_v, n1_h1, n2_v, n2_h1):
    """2D restriction of :func:`_normal_flux_3d` (both components)."""
    theta = rhoY / rho
    f_h1 = mom_h1 * theta
    f_v = mom_v * theta
    return n1_v * f_v + n1_h1 * f_h1, n2_v * f_v + n2_h1 * f_h1


def _momentum_pot_temp_divergence_2d(rho, rhou, rhov, rhoY, dx, dy):
    """2D divergence of the theta-weighted momentum, theta = rhoY/rho."""
    theta = rhoY / rho
    return compute_2d(rhou * theta, rhov * theta, dx, dy)


def _momentum_pot_temp_divergence_3d(rho, rhou, rhov, rhow, rhoY, dx, dy, dz):
    """3D divergence of the theta-weighted momentum, theta = rhoY/rho."""
    theta = rhoY / rho
    return compute_3d_sum(rhou * theta, rhov * theta, rhow * theta, dx, dy, dz)


def wall_zero_dims(ud, ndim):
    """Static tuple of dims whose momentum wall slabs are zeroed (host-side).

    Mirrors the boundary handling at the top of the numpy
    ``compute_at_nodes``: every WALL/RAYLEIGH axis, unless the case uses the
    atmospheric extension.
    """
    if hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        return ()
    return tuple(
        dim
        for dim in range(ndim)
        if ud.bdry_type[dim] in (opts.BdryType.WALL, opts.BdryType.RAYLEIGH)
    )


def compute_at_nodes(momenta, rho, rhoY, ndim, dxyz, wall_dims=(), metric=None):
    """Divergence of the theta-weighted momentum with boundary handling.

    Parameters
    ----------
    momenta : tuple
        (rhou, rhov, rhow); rhow may be None in 2D.
    rho, rhoY : array
    ndim : int (static)
    dxyz : tuple of float
        (dx, dy, dz); dz unused in 2D.
    wall_dims : tuple of int (static)
        Dims whose momentum boundary slabs are zeroed (see
        :func:`wall_zero_dims`).
    metric : None or tuple (static structure)
        Terrain metric ``(N, cart_v, cart_haxes)``: the nested tuple of
        area-normal Cartesian components (outer index = array axis) plus
        the fixed Cartesian role axes — mirrors ``terrain.MetricFields``.

    Returns
    -------
    (rhs, momenta) : the divergence and the slab-zeroed momenta.
    """
    dx, dy, dz = dxyz
    momenta = list(momenta)

    for dim in wall_dims:
        lo, hi = axes.wall_slabs(ndim, dim)
        for i, field in enumerate(momenta):
            if field is None:
                continue
            momenta[i] = field.at[lo].set(0.0).at[hi].set(0.0)
            field = momenta[i]

    if ndim == 2:
        rhou, rhov = momenta[0], momenta[1]
        if metric is not None:
            N, cv, chax = metric
            ch1 = chax[0]
            f_h1, f_v = _normal_fluxes_2d(
                rho, rhoY, rhou, rhov, N[0][cv], N[0][ch1], N[1][cv], N[1][ch1]
            )
            rhs = compute_2d(f_h1, f_v, dx, dy)
        else:
            rhs = _momentum_pot_temp_divergence_2d(rho, rhou, rhov, rhoY, dx, dy)
    elif metric is not None:
        N, cv, chax = metric
        ch1, ch2 = chax
        flux = [
            _normal_flux_3d(
                rho,
                rhoY,
                momenta[cv],
                momenta[ch1],
                momenta[ch2],
                N[a][cv],
                N[a][ch1],
                N[a][ch2],
            )
            for a in range(3)
        ]
        rhs = compute_3d_sum(flux[0], flux[1], flux[2], dx, dy, dz)
    else:
        rhs = _momentum_pot_temp_divergence_3d(
            rho, momenta[0], momenta[1], momenta[2], rhoY, dx, dy, dz
        )

    return rhs, tuple(momenta)
