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


def _metric_contravariant_fluxes(rho, rhoY, mom_h1, mom_v, mom_h2, J, G1, G2):
    """Terrain-following flux components of the theta-weighted momentum."""
    theta = rhoY / rho
    f_h1 = mom_h1 * theta
    f_h2 = mom_h2 * theta
    f_v = mom_v * theta - G1 * f_h1 - G2 * f_h2
    return J * f_h1, f_v, J * f_h2


def _metric_contravariant_fluxes_2d(rho, rhoY, mom_h1, mom_v, J, G1):
    """2D restriction of :func:`_metric_contravariant_fluxes`."""
    theta = rhoY / rho
    f_h1 = mom_h1 * theta
    f_v = mom_v * theta - G1 * f_h1
    return J * f_h1, f_v


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
        Terrain metric (J, G1, G2, vaxis, haxes) with array J/G1/G2;
        G2/haxes ignored in 2D.

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
            J, G1 = metric[0], metric[1]
            f_h1, f_v = _metric_contravariant_fluxes_2d(rho, rhoY, rhou, rhov, J, G1)
            rhs = compute_2d(f_h1, f_v, dx, dy)
        else:
            rhs = _momentum_pot_temp_divergence_2d(rho, rhou, rhov, rhoY, dx, dy)
    elif metric is not None:
        J, G1, G2, vaxis, haxes = metric
        a_h1, a_h2 = haxes
        f_h1, f_v, f_h2 = _metric_contravariant_fluxes(
            rho, rhoY, momenta[a_h1], momenta[vaxis], momenta[a_h2], J, G1, G2
        )
        flux = [None, None, None]
        flux[a_h1], flux[vaxis], flux[a_h2] = f_h1, f_v, f_h2
        rhs = compute_3d_sum(flux[0], flux[1], flux[2], dx, dy, dz)
    else:
        rhs = _momentum_pot_temp_divergence_3d(
            rho, momenta[0], momenta[1], momenta[2], rhoY, dx, dy, dz
        )

    return rhs, tuple(momenta)
