import numpy as np
from numba import njit

from ....utils import options as opts
from ....utils.slices import get_neighbor_indices, get_interface_indices

from ...utils import variable as var

def do(mem, ud, lmbda, split_step, tag=None):
    """
    Reconstruct the limited slopes at the cell interfaces.

    """
    # elem, Sol, flux, th, cache = mem.elem, mem.sol, mem.flux, mem.th, mem.cache
    flux = mem.flux[split_step]
    gamm = mem.th.gamm

    order_two = 1  # always 1

    mem.sol.primitives(mem.th)

    if tag == "rk":
        lmbda = 0.0

    lefts_idx, rights_idx, inner_idx = get_interface_indices(mem.elem.ndim)

    # inner_idx here are where the interface fluxes are calculated with non-zero values.
    face_inner_idx = inner_idx
    u = np.zeros_like(mem.sol.rhoY)
    u[inner_idx] = (
        0.5
        * (flux.rhoY[face_inner_idx][lefts_idx] + flux.rhoY[face_inner_idx][rights_idx])
        / mem.sol.rhoY[inner_idx]
    )

    shape = mem.sol.u.shape

    cache = mem.cache.get_recovery_objects(shape, ud)
    Diffs = cache['Diffs']
    Ampls = cache['Ampls']
    Lefts = cache['Lefts']
    Rights = cache['Rights']

    Diffs.u[..., :-1] = mem.sol.u[rights_idx] - mem.sol.u[lefts_idx]
    Diffs.v[..., :-1] = mem.sol.v[rights_idx] - mem.sol.v[lefts_idx]
    Diffs.w[..., :-1] = mem.sol.w[rights_idx] - mem.sol.w[lefts_idx]
    Diffs.X[..., :-1] = mem.sol.X[rights_idx] - mem.sol.X[lefts_idx]
    Diffs.Y[..., :-1] = 1.0 / mem.sol.Y[rights_idx] - 1.0 / mem.sol.Y[lefts_idx]

    Slopes = slopes(Diffs, ud, mem.elem)

    Ampls.u[...] = 0.5 * Slopes.u * (1.0 - lmbda * u)
    Ampls.v[...] = 0.5 * Slopes.v * (1.0 - lmbda * u)
    Ampls.w[...] = 0.5 * Slopes.w * (1.0 - lmbda * u)
    Ampls.X[...] = 0.5 * Slopes.X * (1.0 - lmbda * u)
    Ampls.Y[...] = 0.5 * Slopes.Y * (1.0 - lmbda * u)

    Lefts.u[...] = mem.sol.u + order_two * Ampls.u
    Lefts.v[...] = mem.sol.v + order_two * Ampls.v
    Lefts.w[...] = mem.sol.w + order_two * Ampls.w
    Lefts.X[...] = mem.sol.X + order_two * Ampls.X
    Lefts.Y[...] = 1.0 / (1.0 / mem.sol.Y + order_two * Ampls.Y)

    Ampls.u[...] = -0.5 * Slopes.u * (1.0 + lmbda * u)
    Ampls.v[...] = -0.5 * Slopes.v * (1.0 + lmbda * u)
    Ampls.w[...] = -0.5 * Slopes.w * (1.0 + lmbda * u)
    Ampls.X[...] = -0.5 * Slopes.X * (1.0 + lmbda * u)
    Ampls.Y[...] = -0.5 * Slopes.Y * (1.0 + lmbda * u)

    Rights.u[...] = mem.sol.u + order_two * Ampls.u
    Rights.v[...] = mem.sol.v + order_two * Ampls.v
    Rights.w[...] = mem.sol.w + order_two * Ampls.w
    Rights.X[...] = mem.sol.X + order_two * Ampls.X
    Rights.Y[...] = 1.0 / (1.0 / mem.sol.Y + order_two * Ampls.Y)

    vel = [mem.sol.u, mem.sol.v, mem.sol.w]

    # Lefts.rhoY[lefts_idx] = Rights.rhoY[rights_idx] = 0.5 * (mem.sol.rhoY[lefts_idx] + mem.sol.rhoY[rights_idx]) \
    #     - order_two * 0.5 * lmbda * (mem.sol.u[rights_idx] * mem.sol.rhoY[rights_idx] - mem.sol.u[lefts_idx] * mem.sol.rhoY[lefts_idx])
    Lefts.rhoY[lefts_idx] = Rights.rhoY[rights_idx] = 0.5 * (
        mem.sol.rhoY[lefts_idx] + mem.sol.rhoY[rights_idx]
    ) - order_two * 0.5 * lmbda * (
        vel[split_step][rights_idx] * mem.sol.rhoY[rights_idx]
        - vel[split_step][lefts_idx] * mem.sol.rhoY[lefts_idx]
    )

    Lefts.p0[lefts_idx] = Rights.p0[rights_idx] = Lefts.rhoY[lefts_idx] ** gamm

    get_conservatives(Rights)
    get_conservatives(Lefts)

    return Lefts, Rights

def slopes(Diffs, ud, elem):
    """Reconstruct piecewise linear slopes in cells."""
    # Configuration
    variable_config = {
        'u': ud.limiter_type_velocity,
        'v': ud.limiter_type_velocity, 
        'w': ud.limiter_type_velocity,
        'X': ud.limiter_type_scalars,
        'Y': ud.limiter_type_scalars
    }
    
    lefts_idx, rights_idx = get_neighbor_indices(elem.ndim)
    
    # Initialize slopes
    # TBD: Move this to initialisation stage
    Slopes = var.Characters(Diffs.u.shape)
    
    # Process each variable
    for var_name, limiter_type in variable_config.items():
        diff_data = getattr(Diffs, var_name)
        
        # Extract amplitudes
        al = diff_data[lefts_idx][lefts_idx]
        ar = diff_data[lefts_idx][rights_idx]
        
        # Calculate and assign slopes
        slope_data = limiters(limiter_type, al, ar)
        getattr(Slopes, var_name)[..., 1:-1] = slope_data
    
    return Slopes

@njit
def limiters(limiter_type, al, ar):
    """
    Applies the limiter type specified in the initial conditions to recovery the slope.

    Parameters
    ----------
    limiter_type : :py:class:`management.enumerator.LimiterType`
        LimiterType list
    al : :py:class:`management.variable.States`
        Left indices of the `Diffs` array for the respective quantities.
    ar : :py:class:`management.variable.States`
        Right indices of the `Diffs` array for the respective quantities.

    Returns
    -------
    :py:class:`management.variable.States`
        The reconstructed slope in the cell

    Attention
    ---------
    For now, only the limiter type `NONE` is supported. This takes $\\frac{(al + ar)}{2}$.
    """
    # write switch for limiter types
    # for now, just use LimiterType == None
    if limiter_type == opts.LimiterType.NONE:
        return 0.5 * (al + ar)


def get_conservatives(U):
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
