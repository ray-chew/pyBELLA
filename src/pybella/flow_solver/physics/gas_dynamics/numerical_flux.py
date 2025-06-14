# -*- coding: utf-8 -*-
import numpy as np

from ....utils.operators import create_convolution_kernels, apply_directional_convolution
from ....utils.slices import get_inner_slice, get_interface_indices, get_last_dim_inner_slice

def recompute_advective_fluxes(mem, **kwargs):
    """Recompute the advective fluxes at the cell interfaces.
    
    Parameters
    ----------
    mem : object
        Memory object containing sol and flux attributes
    **kwargs
        Optional pre-computed velocity components ('u', 'v', 'w')
    """
    ndim = mem.sol.rho.ndim
    inner_idx = get_inner_slice(ndim)
    kernels = create_convolution_kernels(ndim)
    
    # Define the component order and corresponding flux indices
    components = ['u', 'v'] if ndim == 2 else ['u', 'v', 'w']
    rho_components = ['rhou', 'rhov'] if ndim == 2 else ['rhou', 'rhov', 'rhow']
    
    for i, (comp, rho_comp) in enumerate(zip(components, rho_components)):
        # Use provided velocity or compute from momentum
        if comp in kwargs:
            rhoY_vel = kwargs[comp]
        else:
            momentum = getattr(mem.sol, rho_comp)
            rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho
        
        # Apply directional convolution
        mem.flux[i].rhoY[inner_idx] = apply_directional_convolution(
            rhoY_vel, kernels[comp], comp, ndim
        )

def hll_solver(mem, flux, Lefts, Rights):
    """
    HLL solver for the Riemann problem. Chooses the advected quantities from `Lefts` or `Rights` based on the direction given by `flux`.

    Returns
    -------
    :py:class:`management.variable.States`
        `flux` data container with the solution of the Riemann problem.
    
    """
    def _compute_flux_component(flux_attr, state_attr=None, state_value=1.0):
        """Helper function to compute a single flux component."""
        left_weight = upl[left_idx] / Lefts.Y[left_idx]
        right_weight = upr[right_idx] / Rights.Y[right_idx]
        
        if state_attr is not None:
            left_val = getattr(Lefts, state_attr)[left_idx]
            right_val = getattr(Rights, state_attr)[right_idx]
        else:
            left_val = right_val = state_value
        
        getattr(flux, flux_attr)[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (
            left_weight * left_val + right_weight * right_val
        )

    ndim = mem.sol.rho.ndim
    left_idx, right_idx, _ = get_interface_indices(ndim)
    remove_cols_idx = get_last_dim_inner_slice(ndim)

    # Compute primitive variables
    Lefts.primitives(mem.th)
    Rights.primitives(mem.th)

    # Compute upwind weights
    upwind = 0.5 * (1.0 + np.sign(flux.rhoY))
    upl = upwind[right_idx]
    upr = 1.0 - upwind[left_idx]

    # Compute all flux components
    _compute_flux_component('rhou', 'u')
    _compute_flux_component('rho')  # Uses default state_value=1.0
    _compute_flux_component('rhov', 'v')
    _compute_flux_component('rhow', 'w')
    _compute_flux_component('rhoX', 'X')

    return flux