# -*- coding: utf-8 -*-
import numpy as np
from numba import njit

from ....utils.operators import convolution
from ....utils import slices


def recompute(mem, **kwargs):
    """Recompute the advective fluxes at the cell interfaces.

    Parameters
    ----------
    mem : object
        Memory object containing sol and flux attributes
    **kwargs
        Optional pre-computed velocity components ('u', 'v', 'w')
    """
    ndim = mem.sol.rho.ndim
    inner_idx = slices.get_inner_slice(ndim)
    kernels = convolution.get_flux_kernels(ndim)

    # Define the component order and corresponding flux indices
    components = ["u", "v"] if ndim == 2 else ["u", "v", "w"]
    rho_components = ["rhou", "rhov"] if ndim == 2 else ["rhou", "rhov", "rhow"]

    flux = mem.cache.get_flux_containers(mem.elem)

    for i, (comp, rho_comp) in enumerate(zip(components, rho_components)):
        # Use provided velocity or compute from momentum
        if comp in kwargs:
            rhoY_vel = kwargs[comp]
        else:
            momentum = getattr(mem.sol, rho_comp)
            rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho

        # Apply directional convolution
        flux[i].rhoY[inner_idx] = convolution.apply_directional_convolution(
            rhoY_vel, kernels[comp], comp, ndim
        )
