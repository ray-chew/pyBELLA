# -*- coding: utf-8 -*-
import numpy as np
from numba import njit

from ....utils.operators import convolution
from ....utils import slices


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
            if metric is not None and i == metric.vaxis:
                # contravariant vertical mass flux rhoY*(w - G.u_h)/rho
                # (J * eta_dot — what actually crosses an eta-face)
                a_h1, a_h2 = metric.haxes
                momentum = momentum - metric.G1 * getattr(mem.sol, rho_components[a_h1])
                if metric.G2 is not None:
                    momentum = momentum - metric.G2 * getattr(
                        mem.sol, rho_components[a_h2]
                    )
                rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho
            elif metric is not None:
                # horizontal mass fluxes carry the Jacobian (face-area weight)
                rhoY_vel = metric.J * mem.sol.rhoY * momentum / mem.sol.rho
            else:
                rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho

        # Apply directional convolution
        flux[i].rhoY[inner_idx] = convolution.apply_directional_convolution(
            rhoY_vel, kernels[comp], comp, ndim
        )
