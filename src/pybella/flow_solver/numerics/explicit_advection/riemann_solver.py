import numpy as np
import numba as nb
from ....utils import slices


def hll(mem, flux, Lefts, Rights):
    """
    HLL solver for the Riemann problem. Chooses the advected quantities from `Lefts` or `Rights` based on the direction given by `flux`.

    Returns
    -------
    :py:class:`management.variable.States`
        `flux` data container with the solution of the Riemann problem.

    """
    ndim = mem.sol.rho.ndim
    left_idx, right_idx, _ = slices.get_interface_indices(ndim)
    remove_cols_idx = slices.get_last_dim_inner_slice(ndim)

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


@nb.njit(cache=True)
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
