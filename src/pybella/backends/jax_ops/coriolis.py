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
from pybella.flow_solver.numerics import coriolis as coriolis_np


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


@jax.jit
def compute_coefficients_general(wh1, wh2, wv, e1, e2, e3, nu, nonhydro):
    """General (C11) H^-1 coefficients for an arbitrary up-direction e.

    JAX twin of ``numerics.coriolis._compute_coriolis_coefficients_general``
    (Sherman-Morrison on the pure-rotation inverse + the C11 alpha_w
    prefactor). Same FP order as that numba kernel; role symbols
    w = (wh1, wv, wh2), e = (e_h1, e_v, e_h2). Returns (h11..h33, denom).
    """
    w1 = wh1
    w2 = wv
    w3 = wh2
    ooD = 1.0 / (1.0 + (w1 * w1 + w2 * w2 + w3 * w3))

    c11 = (1.0 + w1 * w1) * ooD
    c12 = (w1 * w2 + w3) * ooD
    c13 = (w1 * w3 - w2) * ooD
    c21 = (w2 * w1 - w3) * ooD
    c22 = (1.0 + w2 * w2) * ooD
    c23 = (w2 * w3 + w1) * ooD
    c31 = (w3 * w1 + w2) * ooD
    c32 = (w3 * w2 - w1) * ooD
    c33 = (1.0 + w3 * w3) * ooD

    # t = C^-1 e (column), s = e^T C^-1 (row)
    t1 = c11 * e1 + c12 * e2 + c13 * e3
    t2 = c21 * e1 + c22 * e2 + c23 * e3
    t3 = c31 * e1 + c32 * e2 + c33 * e3
    s1 = e1 * c11 + e2 * c21 + e3 * c31
    s2 = e1 * c12 + e2 * c22 + e3 * c32
    s3 = e1 * c13 + e2 * c23 + e3 * c33

    mu = nu + (nonhydro - 1.0)
    denom = 1.0 / (1.0 + mu * (e1 * t1 + e2 * t2 + e3 * t3))

    h11 = c11 - mu * t1 * s1 * denom
    h12 = c12 - mu * t1 * s2 * denom
    h13 = c13 - mu * t1 * s3 * denom
    h21 = c21 - mu * t2 * s1 * denom
    h22 = c22 - mu * t2 * s2 * denom
    h23 = c23 - mu * t2 * s3 * denom
    h31 = c31 - mu * t3 * s1 * denom
    h32 = c32 - mu * t3 * s2 * denom
    h33 = c33 - mu * t3 * s3 * denom

    # (C11) alpha_w prefactor: G = H^-1 - (1-a)(I - e e^T) H^-1 (e e^T)
    col1 = h11 * e1 + h12 * e2 + h13 * e3
    col2 = h21 * e1 + h22 * e2 + h23 * e3
    col3 = h31 * e1 + h32 * e2 + h33 * e3
    col_par = e1 * col1 + e2 * col2 + e3 * col3
    onema = 1.0 - nonhydro
    cp1 = onema * (col1 - col_par * e1)
    cp2 = onema * (col2 - col_par * e2)
    cp3 = onema * (col3 - col_par * e3)
    h11 = h11 - cp1 * e1
    h12 = h12 - cp1 * e2
    h13 = h13 - cp1 * e3
    h21 = h21 - cp2 * e1
    h22 = h22 - cp2 * e2
    h23 = h23 - cp2 * e3
    h31 = h31 - cp3 * e1
    h32 = h32 - cp3 * e2
    h33 = h33 - cp3 * e3

    return h11, h12, h13, h21, h22, h23, h31, h32, h33, denom


@jax.jit
def apply_inverse_general(U, V, W, wh1, wh2, wv, e1, e2, e3, nu, nonhydro):
    """(u, v, w) = H^-1 @ (U, V, W) with the general up-direction e."""
    h11, h12, h13, h21, h22, h23, h31, h32, h33, _ = compute_coefficients_general(
        wh1, wh2, wv, e1, e2, e3, nu, nonhydro
    )
    u = h11 * U + h12 * V + h13 * W
    v = h21 * U + h22 * V + h23 * W
    w = h31 * U + h32 * V + h33 * W
    return u, v, w


def _role_inputs(mem, ud, dt):
    """Role-ordered (h1, v, h2) H^-1 inputs.

    Reuses the numpy ``role_components`` / ``_up_role_components`` so the
    spatially varying rotation field (``ud.coriolis_field``) and the local
    up-direction ``e`` are built once, cached on ``mem.elem`` and bit-
    identical to the numpy path. ``e_role`` is None on vertical-line/no-
    metric runs (the legacy scalar kernel applies).
    """
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[axes.vertical_axis(ud)]
    Msq = ud.Msq

    ax_h1, ax_v, ax_h2 = axes.role_perm(axes.vertical_axis(ud))
    w_h1, w_v, w_h2 = coriolis_np.role_components(mem, ud)
    wh1, wv, wh2 = dt * w_h1, dt * w_v, dt * w_h2
    strat = mem.npf.HydroState_n.get_dSdy(mem.elem, mem.node)
    Y = mem.sol.rhoY / mem.sol.rho
    nu = -(dt**2) * (g / Msq) * strat * Y
    e_role = coriolis_np._up_role_components(mem, ud)
    return (ax_h1, ax_v, ax_h2), wh1, wv, wh2, nu, float(nonhydro), e_role


def _coeffs(wh1, wh2, wv, nu, nonhydro, e_role):
    """H^-1 coefficient tuple (h11..h33, denom) — general or legacy."""
    if e_role is None:
        return compute_coefficients(
            jnp.asarray(wh1), jnp.asarray(wh2), jnp.asarray(wv), nu, nonhydro
        )
    e1, e2, e3 = (jnp.asarray(e) for e in e_role)
    return compute_coefficients_general(
        jnp.asarray(wh1), jnp.asarray(wh2), jnp.asarray(wv), e1, e2, e3, nu, nonhydro
    )


def multiply_inverse_terms(
    Vec, mem, ud, dt, attrs=("rhou", "rhov", "rhow"), get_coeffs=False
):
    """Drop-in twin of the numpy ``multiply_inverse_terms`` (same contract:
    mutates the ``attrs`` fields of ``Vec`` and optionally returns the 2D
    coefficient block)."""
    (ax_h1, ax_v, ax_h2), wh1, wv, wh2, nu, nonhydro, e_role = _role_inputs(mem, ud, dt)
    nu = jnp.asarray(nu)

    VecU = getattr(Vec, attrs[ax_h1])
    VecV = getattr(Vec, attrs[ax_v])
    VecW = getattr(Vec, attrs[ax_h2])

    if e_role is None:
        u, v, w = apply_inverse(
            jnp.asarray(VecU),
            jnp.asarray(VecV),
            jnp.asarray(VecW),
            jnp.asarray(wh1),
            jnp.asarray(wh2),
            jnp.asarray(wv),
            nu,
            nonhydro,
        )
    else:
        e1, e2, e3 = (jnp.asarray(e) for e in e_role)
        u, v, w = apply_inverse_general(
            jnp.asarray(VecU),
            jnp.asarray(VecV),
            jnp.asarray(VecW),
            jnp.asarray(wh1),
            jnp.asarray(wh2),
            jnp.asarray(wv),
            e1,
            e2,
            e3,
            nu,
            nonhydro,
        )
    VecU[...] = np.asarray(u)
    VecV[...] = np.asarray(v)
    VecW[...] = np.asarray(w)

    if get_coeffs:
        # 2D-only path (the (h1, v) block); fill the shared cache views so
        # downstream consumers see exactly what the numpy twin leaves there
        views = mem.cache.get_coriolis_array_views(nu.shape)
        coeffs = _coeffs(wh1, wh2, wv, nu, nonhydro, e_role)
        for view, val in zip(views, coeffs):
            view[...] = np.asarray(val)
        h11, h12, _, h21, h22 = views[0], views[1], views[2], views[3], views[4]
        return (h11.T, h22.T, h12.T, h21.T)


def compute_inverse_coefficients(mem, ud, dt):
    """Drop-in twin: fill and return the cached role-indexed H^-1 fields."""
    _, wh1, wv, wh2, nu, nonhydro, e_role = _role_inputs(mem, ud, dt)
    views = mem.cache.get_coriolis_array_views(nu.shape)
    coeffs = _coeffs(wh1, wh2, wv, jnp.asarray(nu), nonhydro, e_role)
    for view, val in zip(views, coeffs):
        view[...] = np.asarray(val)
    return views
