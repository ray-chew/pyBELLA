# -*- coding: utf-8 -*-
import numpy as np
from numba import njit

from ....utils.operators import (
    get_flux_convolution_kernels,
    apply_directional_convolution,
)
from ....utils.slices import (
    get_inner_slice,
    get_interface_indices,
    get_last_dim_inner_slice,
)


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
    kernels = get_flux_convolution_kernels(ndim)

    # Define the component order and corresponding flux indices
    components = ["u", "v"] if ndim == 2 else ["u", "v", "w"]
    rho_components = ["rhou", "rhov"] if ndim == 2 else ["rhou", "rhov", "rhow"]

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


@njit(cache=True)
def _compute_flux_component(
    flux_values,
    rhoY_values,
    left_weight,
    right_weight,
    left_val,
    right_val,
    remove_cols_idx,
):
    """Numba-optimised flux component computation."""
    flux_values[remove_cols_idx] = rhoY_values[remove_cols_idx] * (
        left_weight * left_val + right_weight * right_val
    )


def hll_solver(mem, flux, Lefts, Rights):
    """
    HLL solver for the Riemann problem. Chooses the advected quantities from `Lefts` or `Rights` based on the direction given by `flux`.

    Returns
    -------
    :py:class:`management.variable.States`
        `flux` data container with the solution of the Riemann problem.

    """
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

    # Pre-compute common weights
    left_weight = upl[left_idx] / Lefts.Y[left_idx]
    right_weight = upr[right_idx] / Rights.Y[right_idx]

    # Define flux components to compute
    flux_components = [
        ("rhou", "u"),
        ("rho", None, 1.0),  # state_value=1.0
        ("rhov", "v"),
        ("rhow", "w"),
        ("rhoX", "X"),
    ]

    # Compute all flux components
    for component in flux_components:
        flux_attr = component[0]
        state_attr = component[1] if len(component) > 1 else None
        state_value = component[2] if len(component) > 2 else 1.0

        # Get state values
        if state_attr is not None:
            left_val = getattr(Lefts, state_attr)[left_idx]
            right_val = getattr(Rights, state_attr)[right_idx]
        else:
            left_val = right_val = state_value

        _compute_flux_component(
            getattr(flux, flux_attr),
            flux.rhoY,
            left_weight,
            right_weight,
            left_val,
            right_val,
            remove_cols_idx,
        )

    return flux
