import numpy as np

from ....utils import options as opts
from ....utils.slices import get_neighbor_indices

from ...utils import variable as var


def do(mem, ud, lmbda, split_step, tag=None, use_cache=False):
    """
    Reconstruct the limited slopes at the cell interfaces.

    Parameters
    ----------
    Sol : :py:class:`management.variable.Vars`
        Solution data container on cell centers.
    flux : :py:class:`management.variable.States`
        Flux data container on cell interfaces.
    lmbda : float
        :math:`\\frac{dt}{dx}`, where :math:`dx` is the grid-size in the direction of the substep.
    ud : :py:class:`inputs.user_data.UserDataInit`
        Class container for the initial condition.
    th : :py:class:`physics.gas_dynamics.thermodynamic.init`
        Class container for the thermodynamical constants.
    elem : :py:class:`discretization.kgrid.ElemSpaceDiscr`
        Class container for the cell-grid.
    split_step : int
        Tracks the substep in the Strang-splitting.
    tag : `None` or `rk`
        Default is `None` which uses a second-order Strang-splitting. `rk` toggles a first-order Runge-Kutta update for the advection scheme.

    Returns
    -------
    :py:class:`management.variable.States`, :py:class:`management.variable.States`
        Lefts, Rights are containers for the advected quantities at to the left and the right of the cell interfaces.

    """
    elem, Sol, flux, th, cache = mem.elem, mem.sol, mem.flux, mem.th, mem.cache
    flux = flux[split_step]
    gamm = th.gamm

    order_two = 1  # always 1

    Sol.primitives(th)

    if tag == "rk":
        lmbda = 0.0

    ndim = elem.ndim
    lefts_idx, rights_idx, inner_idx = (
        [
            slice(
                None,
            )
        ]
        * ndim,
        [
            slice(
                None,
            )
        ]
        * ndim,
        [slice(1, -1)] * ndim,
    )
    lefts_idx[-1] = slice(0, -1)
    rights_idx[-1] = slice(1, None)
    lefts_idx, rights_idx, inner_idx = (
        tuple(lefts_idx),
        tuple(rights_idx),
        tuple(inner_idx),
    )

    # inner_idx here are where the interface fluxes are calculated with non-zero values.
    face_inner_idx = inner_idx
    u = np.zeros_like(Sol.rhoY)
    u[inner_idx] = (
        0.5
        * (flux.rhoY[face_inner_idx][lefts_idx] + flux.rhoY[face_inner_idx][rights_idx])
        / Sol.rhoY[inner_idx]
    )

    shape = Sol.u.shape
    
    if use_cache:
        cache = cache.get_recovery_objects(shape, ud)
        Diffs = cache['Diffs']
        Ampls = cache['Ampls']
        Lefts = cache['Lefts']
        Rights = cache['Rights']
    else:
        # Fallback
        Diffs = var.States(shape, ud)
        Ampls = var.Characters(shape)
        Lefts = var.States(shape, ud)
        Rights = var.States(shape, ud)

    Diffs.u[..., :-1] = Sol.u[rights_idx] - Sol.u[lefts_idx]
    Diffs.v[..., :-1] = Sol.v[rights_idx] - Sol.v[lefts_idx]
    Diffs.w[..., :-1] = Sol.w[rights_idx] - Sol.w[lefts_idx]
    Diffs.X[..., :-1] = Sol.X[rights_idx] - Sol.X[lefts_idx]
    Diffs.Y[..., :-1] = 1.0 / Sol.Y[rights_idx] - 1.0 / Sol.Y[lefts_idx]

    Slopes = slopes(Diffs, ud, elem)

    Ampls.u[...] = 0.5 * Slopes.u * (1.0 - lmbda * u)
    Ampls.v[...] = 0.5 * Slopes.v * (1.0 - lmbda * u)
    Ampls.w[...] = 0.5 * Slopes.w * (1.0 - lmbda * u)
    Ampls.X[...] = 0.5 * Slopes.X * (1.0 - lmbda * u)
    Ampls.Y[...] = 0.5 * Slopes.Y * (1.0 - lmbda * u)

    Lefts.u[...] = Sol.u + order_two * Ampls.u
    Lefts.v[...] = Sol.v + order_two * Ampls.v
    Lefts.w[...] = Sol.w + order_two * Ampls.w
    Lefts.X[...] = Sol.X + order_two * Ampls.X
    Lefts.Y[...] = 1.0 / (1.0 / Sol.Y + order_two * Ampls.Y)

    Ampls.u[...] = -0.5 * Slopes.u * (1.0 + lmbda * u)
    Ampls.v[...] = -0.5 * Slopes.v * (1.0 + lmbda * u)
    Ampls.w[...] = -0.5 * Slopes.w * (1.0 + lmbda * u)
    Ampls.X[...] = -0.5 * Slopes.X * (1.0 + lmbda * u)
    Ampls.Y[...] = -0.5 * Slopes.Y * (1.0 + lmbda * u)

    Rights.u[...] = Sol.u + order_two * Ampls.u
    Rights.v[...] = Sol.v + order_two * Ampls.v
    Rights.w[...] = Sol.w + order_two * Ampls.w
    Rights.X[...] = Sol.X + order_two * Ampls.X
    Rights.Y[...] = 1.0 / (1.0 / Sol.Y + order_two * Ampls.Y)

    vel = [Sol.u, Sol.v, Sol.w]

    # Lefts.rhoY[lefts_idx] = Rights.rhoY[rights_idx] = 0.5 * (Sol.rhoY[lefts_idx] + Sol.rhoY[rights_idx]) \
    #     - order_two * 0.5 * lmbda * (Sol.u[rights_idx] * Sol.rhoY[rights_idx] - Sol.u[lefts_idx] * Sol.rhoY[lefts_idx])
    Lefts.rhoY[lefts_idx] = Rights.rhoY[rights_idx] = 0.5 * (
        Sol.rhoY[lefts_idx] + Sol.rhoY[rights_idx]
    ) - order_two * 0.5 * lmbda * (
        vel[split_step][rights_idx] * Sol.rhoY[rights_idx]
        - vel[split_step][lefts_idx] * Sol.rhoY[lefts_idx]
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
