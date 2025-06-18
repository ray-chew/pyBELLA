import numpy as np
from numba import njit

from ....utils import options as opts
from ....utils.slices import get_neighbor_indices, get_interface_indices


def compute(mem, flux, ud, lmbda, split_step, tag=None):
    """
    Reconstruct the limited slopes at the cell interfaces.

    """
    gamm = mem.th.gamm

    order_two = 1  # always 1

    mem.sol.primitives(mem.th)

    if tag == "rk":
        lmbda = 0.0

    # TBD: Consider moving these to the cache
    lefts_idx, rights_idx, face_inner_idx = get_interface_indices(mem.elem.ndim)

    # inner_idx here are where the interface fluxes are calculated with non-zero values.
    # TBD: Move to the cache
    u = np.zeros_like(mem.sol.rhoY)
    u[face_inner_idx] = (
        0.5
        * (flux.rhoY[face_inner_idx][lefts_idx] + flux.rhoY[face_inner_idx][rights_idx])
        / mem.sol.rhoY[face_inner_idx]
    )

    shape = mem.sol.u.shape

    # Get cached objects
    cache = mem.cache.get_recovery_objects(shape, ud)
    Diffs, Ampls, Lefts, Rights, Slopes = (
        cache["Diffs"],
        cache["Ampls"],
        cache["Lefts"],
        cache["Rights"],
        cache["Slopes"],
    )

    # Compute differences
    _compute_differences(mem.sol, rights_idx, lefts_idx, Diffs)

    # Compute slopes
    slopes_obj = _slopes(Diffs, Slopes, ud, mem.elem)

    # Compute left-side amplitudes and values
    _compute_amplitudes(
        slopes_obj, lmbda, u, Ampls, sign_factor=1.0, lambda_factor=-1.0
    )
    _compute_reconstructed_values(mem.sol, Ampls, order_two, Lefts)

    # Compute right-side amplitudes and values
    _compute_amplitudes(
        slopes_obj, lmbda, u, Ampls, sign_factor=-1.0, lambda_factor=1.0
    )
    _compute_reconstructed_values(mem.sol, Ampls, order_two, Rights)

    # Return velocity components
    vel = [mem.sol.u, mem.sol.v, mem.sol.w]

    # Compute rhoY reconstruction
    reconstructed_rhoy = _compute_rhoy_reconstruction(
        mem, Lefts, Rights, lefts_idx, rights_idx, vel, split_step, order_two, lmbda
    )

    # Compute pressure reconstruction
    _compute_pressure_reconstruction(
        Lefts, Rights, lefts_idx, rights_idx, reconstructed_rhoy, gamm
    )

    _get_conservatives(Rights)
    _get_conservatives(Lefts)

    return Lefts, Rights


def _slopes(Diffs, Slopes, ud, elem):
    """Reconstruct piecewise linear slopes in cells."""
    # Configuration
    variable_config = {
        "u": ud.limiter_type_velocity,
        "v": ud.limiter_type_velocity,
        "w": ud.limiter_type_velocity,
        "X": ud.limiter_type_scalars,
        "Y": ud.limiter_type_scalars,
    }

    # TBD: Consider moving indices to cache
    lefts_idx, rights_idx = get_neighbor_indices(elem.ndim)

    # Process each variable
    for var_name, limiter_type in variable_config.items():
        diff_data = getattr(Diffs, var_name)

        # Extract amplitudes
        al = diff_data[lefts_idx][lefts_idx]
        ar = diff_data[lefts_idx][rights_idx]

        # Calculate and assign slopes
        slope_data = _limiters(limiter_type, al, ar)
        getattr(Slopes, var_name)[..., 1:-1] = slope_data

    return Slopes


@njit(cache=True)
def _limiters(limiter_type, al, ar):
    """
    Applies the limiter type specified in the initial conditions to recovery the slope.

    """
    # write switch for limiter types
    # for now, just use LimiterType == None
    if limiter_type == opts.LimiterType.NONE:
        return 0.5 * (al + ar)


def _get_conservatives(U):
    """
    Get advected (conservative) quantities at the left and right of the cell interfaces.

    Parameters
    ----------
    U : :py:class:`management.variable.States`
        `Lefts` and `Rights` corresponding to the values at the cell interfaces.
    """
    U.rho = U.rhoY / U.Y
    U.rhou = U.u * U.rho
    U.rhov = U.v * U.rho
    U.rhow = U.w * U.rho
    U.rhoY = U.Y * U.rho
    U.rhoX = U.X * U.rho


def _compute_differences(sol, rights_idx, lefts_idx, diffs):
    """Compute differences between right and left indices for all fields."""
    fields = ["u", "v", "w", "X"]
    for field in fields:
        getattr(diffs, field)[..., :-1] = (
            getattr(sol, field)[rights_idx] - getattr(sol, field)[lefts_idx]
        )

    # Y field has special handling (reciprocal differences)
    diffs.Y[..., :-1] = 1.0 / sol.Y[rights_idx] - 1.0 / sol.Y[lefts_idx]


def _compute_amplitudes(slopes, lmbda, u, ampls, sign_factor=1.0, lambda_factor=1.0):
    """Compute amplitudes for all fields with given sign and lambda factors."""
    fields = ["u", "v", "w", "X", "Y"]
    factor = sign_factor * 0.5 * (1.0 + lambda_factor * lmbda * u)

    for field in fields:
        getattr(ampls, field)[...] = factor * getattr(slopes, field)


def _compute_reconstructed_values(sol, ampls, order_two, result):
    """Compute reconstructed values for all fields."""
    fields = ["u", "v", "w", "X"]
    for field in fields:
        getattr(result, field)[...] = getattr(sol, field) + order_two * getattr(
            ampls, field
        )

    # Y field has special handling (reciprocal computation)
    result.Y[...] = 1.0 / (1.0 / sol.Y + order_two * ampls.Y)


def _compute_rhoy_reconstruction(
    mem, lefts, rights, lefts_idx, rights_idx, vel, split_step, order_two, lmbda
):
    """Compute rhoY reconstruction for both left and right sides."""
    # Use existing arrays for in-place operations
    # First compute on lefts.rhoY, then copy to rights.rhoY

    # Start with average in lefts array
    lefts.rhoY[lefts_idx] = 0.5 * (mem.sol.rhoY[lefts_idx] + mem.sol.rhoY[rights_idx])

    # Subtract velocity correction in-place
    lefts.rhoY[lefts_idx] -= (
        order_two
        * 0.5
        * lmbda
        * (
            vel[split_step][rights_idx] * mem.sol.rhoY[rights_idx]
            - vel[split_step][lefts_idx] * mem.sol.rhoY[lefts_idx]
        )
    )

    # Copy result to rights
    rights.rhoY[rights_idx] = lefts.rhoY[lefts_idx]

    return lefts.rhoY[lefts_idx]


def _compute_pressure_reconstruction(
    lefts, rights, lefts_idx, rights_idx, reconstructed_rhoy, gamm
):
    """Compute pressure reconstruction using power law."""
    reconstructed_p0 = reconstructed_rhoy**gamm
    lefts.p[lefts_idx] = reconstructed_p0
    rights.p[rights_idx] = reconstructed_p0

    return reconstructed_p0
