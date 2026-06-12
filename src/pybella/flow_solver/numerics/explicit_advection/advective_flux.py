# -*- coding: utf-8 -*-
import numpy as np
from numba import njit

from ....utils.operators import convolution
from ....utils import slices

_MOMENTA = ("rhou", "rhov", "rhow")


def _normal_momentum(sol, metric, i):
    """Contravariant face momentum N_i . m for sweep axis i.

    Cartesian-component contraction, vertical component first — the
    bit-exact reduction order of the Phase-0 contract. Shared by the
    numpy and JAX flux assemblies (plain elementwise numpy either way).
    """
    Ni = metric.N[i]
    cv = metric.cart_v
    ch1, ch2 = metric.cart_haxes
    out = Ni[cv] * getattr(sol, _MOMENTA[cv]) + Ni[ch1] * getattr(sol, _MOMENTA[ch1])
    if ch2 is not None:
        out = out + Ni[ch2] * getattr(sol, _MOMENTA[ch2])
    return out


def recompute(mem, ud=None, **kwargs):
    """Recompute the advective fluxes at the cell interfaces.

    Parameters
    ----------
    mem : object
        Memory object containing sol and flux attributes
    ud : UserDataInit, optional
        When given and ud.backend == "jax", the directional convolution
        runs on the JAX backend.
    **kwargs
        Optional pre-computed velocity components ('u', 'v', 'w')
    """
    if ud is not None and getattr(ud, "backend", "numpy") in ("jax", "jax-device"):
        from ....backends.jax_ops import advection as jax_advection

        return jax_advection.recompute_advective_flux(mem, **kwargs)

    ndim = mem.sol.rho.ndim
    inner_idx = slices.get_inner_slice(ndim)
    kernels = convolution.get_flux_kernels(ndim)

    # Define the component order and corresponding flux indices
    components = ["u", "v"] if ndim == 2 else ["u", "v", "w"]
    rho_components = ["rhou", "rhov"] if ndim == 2 else ["rhou", "rhov", "rhow"]

    flux = mem.cache.get_flux_containers(mem.elem)

    # terrain metric; recompute is only called in the unflipped orientation
    # (time_update, between sweeps), so component i matches array axis i
    metric = mem.elem.metric

    for i, (comp, rho_comp) in enumerate(zip(components, rho_components)):
        # Use provided velocity or compute from momentum
        if comp in kwargs:
            rhoY_vel = kwargs[comp]
        else:
            momentum = getattr(mem.sol, rho_comp)
            if metric is not None:
                # general curvilinear mass flux rhoY * (N_i . m) / rho
                # (J * xi_i-dot * rhoY — what actually crosses a xi_i-face);
                # vertical-first contraction so the vertical sweep reduces
                # bit-exactly to the legacy contravariant flux and the
                # horizontals to the J-weighted one
                momentum = _normal_momentum(mem.sol, metric, i)
                rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho
            else:
                rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho

        # Apply directional convolution
        flux[i].rhoY[inner_idx] = convolution.apply_directional_convolution(
            rhoY_vel, kernels[comp], comp, ndim
        )
