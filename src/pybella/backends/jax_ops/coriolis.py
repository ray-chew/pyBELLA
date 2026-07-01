"""JAX twin of :mod:`pybella.flow_solver.numerics.coriolis`.

The two numba kernels become pure jitted functions; the mem-level wrappers
keep the numpy twins' signatures and write results back into the same
containers (Vec attributes, coriolis cache views), so they are drop-in at
the backend seam. Role/axis plumbing is reused from ``utils.axes``.
"""

import functools

import numpy as np
import jax
import jax.numpy as jnp

from pybella.utils import axes


@jax.jit
def compute_coefficients(wh1, wh2, wv, nu, nonhydro):
    """H^-1 coefficients (eq. C11); returns (h11..h33, denom)."""
    wh1_sq = wh1 * wh1
    wh2_sq = wh2 * wh2
    wv_sq = wv * wv
    nu_nh = nu + nonhydro

    denom = 1.0 / (wh1_sq + wh2_sq + nu_nh * (wv_sq + 1.0))

    h11 = (wh1_sq + nu_nh) * denom
    h12 = nonhydro * (wh1 * wv + wh2) * denom
    h13 = (wh1 * wh2 - nu_nh * wv) * denom

    h21 = (wh1 * wv - wh2) * denom
    # h22 carries NO nonhydro factor: see numpy twin
    # (numerics/coriolis._compute_coriolis_coefficients) and
    # dev_notes/hydrostatic_blending.md, Phase H1a. Bit-identical for alpha_w=1.
    h22 = (1.0 + wv_sq) * denom
    h23 = (wh2 * wv + wh1) * denom

    h31 = (wh1 * wh2 + nu_nh * wv) * denom
    h32 = nonhydro * (wh2 * wv - wh1) * denom
    h33 = (nu_nh + wh2_sq) * denom

    return h11, h12, h13, h21, h22, h23, h31, h32, h33, denom


@jax.jit
def apply_inverse(U, V, W, wh1, wh2, wv, nu, nonhydro):
    """(u, v, w) = H^-1 @ (U, V, W), role-ordered (h1, v, h2)."""
    h11, h12, h13, h21, h22, h23, h31, h32, h33, _ = compute_coefficients(
        wh1, wh2, wv, nu, nonhydro
    )
    u = h11 * U + h12 * V + h13 * W
    v = h21 * U + h22 * V + h23 * W
    w = h31 * U + h32 * V + h33 * W
    return u, v, w


def _role_inputs(mem, ud, dt):
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[axes.vertical_axis(ud)]
    Msq = ud.Msq

    ax_h1, ax_v, ax_h2 = axes.role_perm(axes.vertical_axis(ud))
    wdt = dt * ud.coriolis_strength
    wh1, wv, wh2 = wdt[ax_h1], wdt[ax_v], wdt[ax_h2]
    strat = mem.npf.HydroState_n.get_dSdy(mem.elem, mem.node)
    Y = mem.sol.rhoY / mem.sol.rho
    nu = -(dt**2) * (g / Msq) * strat * Y
    return (ax_h1, ax_v, ax_h2), wh1, wv, wh2, nu, float(nonhydro)


def multiply_inverse_terms(
    Vec, mem, ud, dt, attrs=("rhou", "rhov", "rhow"), get_coeffs=False
):
    """Drop-in twin of the numpy ``multiply_inverse_terms`` (same contract:
    mutates the ``attrs`` fields of ``Vec`` and optionally returns the 2D
    coefficient block)."""
    (ax_h1, ax_v, ax_h2), wh1, wv, wh2, nu, nonhydro = _role_inputs(mem, ud, dt)

    VecU = getattr(Vec, attrs[ax_h1])
    VecV = getattr(Vec, attrs[ax_v])
    VecW = getattr(Vec, attrs[ax_h2])

    u, v, w = apply_inverse(
        jnp.asarray(VecU),
        jnp.asarray(VecV),
        jnp.asarray(VecW),
        wh1,
        wh2,
        wv,
        jnp.asarray(nu),
        nonhydro,
    )
    VecU[...] = np.asarray(u)
    VecV[...] = np.asarray(v)
    VecW[...] = np.asarray(w)

    if get_coeffs:
        # 2D-only path (the (h1, v) block); fill the shared cache views so
        # downstream consumers see exactly what the numpy twin leaves there
        views = mem.cache.get_coriolis_array_views(nu.shape)
        coeffs = compute_coefficients(wh1, wh2, wv, jnp.asarray(nu), nonhydro)
        for view, val in zip(views, coeffs):
            view[...] = np.asarray(val)
        h11, h12, _, h21, h22 = views[0], views[1], views[2], views[3], views[4]
        return (h11.T, h22.T, h12.T, h21.T)


def compute_inverse_coefficients(mem, ud, dt):
    """Drop-in twin: fill and return the cached role-indexed H^-1 fields."""
    _, wh1, wv, wh2, nu, nonhydro = _role_inputs(mem, ud, dt)
    views = mem.cache.get_coriolis_array_views(nu.shape)
    coeffs = compute_coefficients(wh1, wh2, wv, jnp.asarray(nu), nonhydro)
    for view, val in zip(views, coeffs):
        view[...] = np.asarray(val)
    return views
