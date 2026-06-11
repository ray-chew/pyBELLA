"""JAX twin of the explicit-advection flux computation (recovery + HLL).

Seam design: the numpy driver in ``explicit_advection.compute_advection``
(dimensional sweeps, array flips, ghost-cell fills, flux-difference
updates) is shared between backends — only the per-sweep flux computation
(`recovery.compute` + `riemann_solver.hll`, the hot kernel) is swapped for
the jitted :func:`_recovery_hll` here. The persistent flux containers are
updated with the same partial writes as the numpy path (HLL touches only
the inner faces; the advective rhoY flux only the inner box), so container
state across sweeps is bit-faithful.

`mem.sol.primitives` is invoked host-side (it is plain numpy) exactly as
the numpy path does, so its side effect on ``sol.u/v/w/Y/X/p`` — and the
values the kernel consumes — are identical by construction.

Limiters: only ``LimiterType.NONE`` exists in the numpy path (central
average); asserted here.
"""

import functools

import numpy as np
import jax
import jax.numpy as jnp

from pybella.utils import options as opts
from pybella.utils import slices
from pybella.utils.operators.convolution import get_flux_kernels

from . import convolution as jax_convolution


@functools.partial(jax.jit, static_argnames=("ndim", "use_metric"))
def _recovery_hll(
    u,
    v,
    w,
    X,
    Y,
    rhoY,
    vel_split,
    flux_rhoY,
    flux_rho_old,
    flux_rhou_old,
    flux_rhov_old,
    flux_rhow_old,
    flux_rhoX_old,
    ooJ,
    lmbda,
    gamm,
    ndim,
    use_metric,
):
    lefts_idx, rights_idx, face_inner_idx = slices.get_interface_indices(ndim)
    remove_cols_idx = slices.get_last_dim_inner_slice(ndim)

    # ---- interface Courant velocity (recovery.compute lines 27-37) ----
    u_face = jnp.zeros_like(rhoY)
    u_face = u_face.at[face_inner_idx].set(
        0.5
        * (flux_rhoY[face_inner_idx][lefts_idx] + flux_rhoY[face_inner_idx][rights_idx])
        / rhoY[face_inner_idx]
    )
    if use_metric:
        u_face = u_face * ooJ

    # ---- differences and central (NONE-limiter) slopes ----
    prims = {"u": u, "v": v, "w": w, "X": X}
    slopes = {}
    for name, field in prims.items():
        d = field[rights_idx] - field[lefts_idx]
        slopes[name] = (
            jnp.zeros_like(field).at[..., 1:-1].set(0.5 * (d[..., :-1] + d[..., 1:]))
        )
    dY = 1.0 / Y[rights_idx] - 1.0 / Y[lefts_idx]
    slopes["Y"] = (
        jnp.zeros_like(Y).at[..., 1:-1].set(0.5 * (dY[..., :-1] + dY[..., 1:]))
    )

    # ---- left/right reconstruction (amplitudes, order_two = 1) ----
    def reconstruct(sign_factor, lambda_factor):
        factor = sign_factor * 0.5 * (1.0 + lambda_factor * lmbda * u_face)
        out = {name: prims[name] + factor * slopes[name] for name in prims}
        out["Y"] = 1.0 / (1.0 / Y + factor * slopes["Y"])
        return out

    Lefts = reconstruct(1.0, -1.0)
    Rights = reconstruct(-1.0, 1.0)

    # ---- rhoY reconstruction (recovery lines 187-201): zeros at the
    # untouched tail like a fresh container; the tail is never read ----
    rhoy_rec = 0.5 * (rhoY[lefts_idx] + rhoY[rights_idx]) - 0.5 * lmbda * (
        vel_split[rights_idx] * rhoY[rights_idx]
        - vel_split[lefts_idx] * rhoY[lefts_idx]
    )
    Lefts["rhoY"] = jnp.zeros_like(rhoY).at[lefts_idx].set(rhoy_rec)
    Rights["rhoY"] = jnp.zeros_like(rhoY).at[rights_idx].set(rhoy_rec)

    # ---- conservatives at the interfaces (recovery _get_conservatives) ----
    def conservatives(U):
        rho = U["rhoY"] / U["Y"]
        return {
            "rho": rho,
            "rhou": U["u"] * rho,
            "rhov": U["v"] * rho,
            "rhow": U["w"] * rho,
            "rhoY": U["Y"] * rho,
            "rhoX": U["X"] * rho,
        }

    Lc = conservatives(Lefts)
    Rc = conservatives(Rights)

    # ---- HLL upwind fluxes (riemann_solver.hll) ----
    def primitives(C):
        return {
            "u": C["rhou"] / C["rho"],
            "v": C["rhov"] / C["rho"],
            "w": C["rhow"] / C["rho"],
            "Y": C["rhoY"] / C["rho"],
            "X": C["rhoX"] / C["rho"],
        }

    Lp = primitives(Lc)
    Rp = primitives(Rc)

    upwind = 0.5 * (1.0 + jnp.sign(flux_rhoY))
    upl = upwind[rights_idx]
    upr = 1.0 - upwind[lefts_idx]

    left_weight = upl[lefts_idx] / Lp["Y"][lefts_idx]
    right_weight = upr[rights_idx] / Rp["Y"][rights_idx]

    def flux_component(old, left_val, right_val):
        return old.at[remove_cols_idx].set(
            flux_rhoY[remove_cols_idx]
            * (left_weight * left_val + right_weight * right_val)
        )

    flux_rhou = flux_component(flux_rhou_old, Lp["u"][lefts_idx], Rp["u"][rights_idx])
    flux_rho = flux_component(flux_rho_old, 1.0, 1.0)
    flux_rhov = flux_component(flux_rhov_old, Lp["v"][lefts_idx], Rp["v"][rights_idx])
    flux_rhow = flux_component(flux_rhow_old, Lp["w"][lefts_idx], Rp["w"][rights_idx])
    flux_rhoX = flux_component(flux_rhoX_old, Lp["X"][lefts_idx], Rp["X"][rights_idx])

    return flux_rho, flux_rhou, flux_rhov, flux_rhow, flux_rhoX


def compute_flux(mem, flux, ud, lmbda, split_step, tag=None):
    """Drop-in twin of recovery.compute + riemann_solver.hll for one sweep.

    Reads the sweep-oriented ``mem.sol`` and ``flux.rhoY``, fills the flux
    container's component fields in place (inner faces only, like the numpy
    path), and returns the container.
    """
    assert ud.limiter_type_velocity == opts.LimiterType.NONE
    assert ud.limiter_type_scalars == opts.LimiterType.NONE

    # same host-side side effect as the numpy recovery (plain numpy)
    mem.sol.primitives(mem.th)

    lmbda_rec = 0.0 if tag == "rk" else lmbda
    use_metric = mem.elem.metric is not None
    ooJ = mem.elem.metric.ooJ if use_metric else np.float64(1.0)
    vel = [mem.sol.u, mem.sol.v, mem.sol.w][split_step]

    out = _recovery_hll(
        jnp.asarray(mem.sol.u),
        jnp.asarray(mem.sol.v),
        jnp.asarray(mem.sol.w),
        jnp.asarray(mem.sol.X),
        jnp.asarray(mem.sol.Y),
        jnp.asarray(mem.sol.rhoY),
        jnp.asarray(vel),
        jnp.asarray(flux.rhoY),
        jnp.asarray(flux.rho),
        jnp.asarray(flux.rhou),
        jnp.asarray(flux.rhov),
        jnp.asarray(flux.rhow),
        jnp.asarray(flux.rhoX),
        jnp.asarray(ooJ),
        lmbda_rec,
        mem.th.gamm,
        mem.elem.ndim,
        use_metric,
    )
    flux_rho, flux_rhou, flux_rhov, flux_rhow, flux_rhoX = out
    flux.rho[...] = np.asarray(flux_rho)
    flux.rhou[...] = np.asarray(flux_rhou)
    flux.rhov[...] = np.asarray(flux_rhov)
    flux.rhow[...] = np.asarray(flux_rhow)
    flux.rhoX[...] = np.asarray(flux_rhoX)

    return flux


def recompute_advective_flux(mem, **kwargs):
    """Drop-in twin of ``advective_flux.recompute``: the rhoY_vel assembly
    stays host-side numpy (bit-identical to the numpy path); the directional
    convolution runs on the JAX backend."""
    ndim = mem.sol.rho.ndim
    inner_idx = slices.get_inner_slice(ndim)
    kernels = get_flux_kernels(ndim)

    components = ["u", "v"] if ndim == 2 else ["u", "v", "w"]
    rho_components = ["rhou", "rhov"] if ndim == 2 else ["rhou", "rhov", "rhow"]

    flux = mem.cache.get_flux_containers(mem.elem)
    metric = mem.elem.metric

    for i, (comp, rho_comp) in enumerate(zip(components, rho_components)):
        if comp in kwargs:
            rhoY_vel = kwargs[comp]
        else:
            momentum = getattr(mem.sol, rho_comp)
            if metric is not None and i == metric.vaxis:
                a_h1, a_h2 = metric.haxes
                momentum = momentum - metric.G1 * getattr(mem.sol, rho_components[a_h1])
                if metric.G2 is not None:
                    momentum = momentum - metric.G2 * getattr(
                        mem.sol, rho_components[a_h2]
                    )
                rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho
            elif metric is not None:
                rhoY_vel = metric.J * mem.sol.rhoY * momentum / mem.sol.rho
            else:
                rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho

        flux[i].rhoY[inner_idx] = np.asarray(
            jax_convolution.apply_directional_convolution(
                rhoY_vel, kernels[comp], comp, ndim
            )
        )
