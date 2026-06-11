"""JAX twin of :func:`...boundary.rayleigh_boundary.rayleigh_damping`.

The sponge profiles (ud.tcy/tny, ud.forcing_tcy/tny) are precomputed 1D
arrays (or 0.0 sentinels); the kernel is pure elementwise relaxation.
``apply_rayleigh_forcing`` stays a host orchestrator on the numpy path —
it evaluates the forcing eigenfunction host-side (t is host-known) and
calls ``rayleigh_damping`` below through the seam.
"""

import functools

import numpy as np
import jax
import jax.numpy as jnp

from pybella.utils import axes
from pybella.utils import options as opts


@functools.partial(jax.jit, static_argnames=("ndim", "has_forcing"))
def _damp(
    rho,
    rhou,
    rhov,
    rhow,
    rhoY,
    p2_nodes,
    tcy,
    tcy_f,
    tny_f,
    Ybar,
    u_wind,
    v_wind,
    w_wind,
    u_f,
    v_f,
    Y_f,
    pi_f,
    mfac,
    c_f,
    ndim,
    has_forcing,
):
    u = rhou / rho
    v = rhov / rho
    Y = rhoY / rho

    if has_forcing:
        p2_nodes = p2_nodes + (tny_f * p2_nodes + jnp.abs(tny_f) * mfac * pi_f)

    u = u + (
        tcy * (u - u_wind) + c_f * (tcy_f * (u - u_wind) + jnp.abs(tcy_f) * mfac * u_f)
    )
    v = v + (
        tcy * (v - v_wind) + c_f * (tcy_f * (v - v_wind) + jnp.abs(tcy_f) * mfac * v_f)
    )
    Y = Y + (
        tcy * (Y - Ybar) + c_f * (tcy_f * (Y - Ybar) + jnp.abs(tcy_f) * mfac * Y_f)
    )

    rhou = rho * u
    rhov = rho * v
    rhoY = rho * Y

    if ndim == 3:
        w = rhow / rho
        w = w + (tcy * (w - w_wind) + c_f * tcy_f * (w - w_wind))
        rhow = rho * w

    return rhou, rhov, rhow, rhoY, p2_nodes


def _vertical_profile(profile, ndim, vaxis):
    if not isinstance(profile, np.ndarray) or profile.ndim != 1:
        return profile
    shape = [1] * ndim
    shape[vaxis] = -1
    return profile.reshape(shape)


def rayleigh_damping(sol, npf, ud, forcing=None):
    """Drop-in twin of the numpy rayleigh_damping (mutates sol, npf)."""
    ndim = sol.rho.ndim
    vaxis = axes.vertical_axis(ud)

    if ud.bdry_type[vaxis] == opts.BdryType.RAYLEIGH:
        tcy = ud.tcy
    else:
        tcy = 0.0
    tcy = _vertical_profile(tcy, ndim, vaxis)

    if forcing is not None:
        tcy_f = _vertical_profile(ud.forcing_tcy, ndim, vaxis)
        tny_f = ud.forcing_tny
        tcy = 0.0
        u_f, v_f, Y_f, pi_f, t = forcing
        if ud.rayleigh_forcing_type == "file":
            raise NotImplementedError(
                "file-based rayleigh forcing is not supported on the JAX "
                "backend; use the numpy backend"
            )
        mfac = 1.0
        c_f = 1.0
    else:
        u_f, v_f, Y_f, pi_f = 0.0, 0.0, 0.0, 0.0
        tcy_f, tny_f = 0.0, 0.0
        mfac = 0.0
        c_f = 0.0

    if npf.HydroState.field_mode:
        Ybar = npf.HydroState.Y0
    else:
        Ybar = _vertical_profile(npf.HydroState.Y0, ndim, vaxis)

    out = _damp(
        jnp.asarray(sol.rho),
        jnp.asarray(sol.rhou),
        jnp.asarray(sol.rhov),
        jnp.asarray(sol.rhow),
        jnp.asarray(sol.rhoY),
        jnp.asarray(npf.p2_nodes),
        jnp.asarray(tcy),
        jnp.asarray(tcy_f),
        jnp.asarray(tny_f),
        jnp.asarray(Ybar),
        ud.u_wind_speed,
        ud.v_wind_speed,
        ud.w_wind_speed,
        jnp.asarray(u_f),
        jnp.asarray(v_f),
        jnp.asarray(Y_f),
        jnp.asarray(pi_f),
        mfac,
        c_f,
        ndim,
        forcing is not None,
    )
    rhou, rhov, rhow, rhoY, p2_nodes = out
    sol.rhou[...] = np.asarray(rhou)
    sol.rhov[...] = np.asarray(rhov)
    sol.rhoY[...] = np.asarray(rhoY)
    if ndim == 3:
        sol.rhow[...] = np.asarray(rhow)
    if forcing is not None:
        npf.p2_nodes[...] = np.asarray(p2_nodes)
