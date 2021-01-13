import numpy as np

from management.variable import States, Characters
from management.enumerator import LimiterType
from .eos import rhoe

def recovery(Sol, flux, lmbda, ud, th, elem, split_step, tag):
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
    th : :py:class:`physics.gas_dynamics.thermodynamic.ThermodynamicInit`
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
    gamm = th.gamm
    
    order_two = 1 # always 1

    Sol.primitives(th)

    # lefts_idx = (slice(None),slice(0,-1))
    # rights_idx = (slice(None),slice(1,None))
    # inner_idx = (slice(1,-1), slice(1,-1))

    if tag == 'rk':
        lmbda = 0.0

    ndim = elem.ndim
    lefts_idx, rights_idx, inner_idx = [slice(None,)] * ndim, [slice(None,)] * ndim, [slice(1,-1)] * ndim
    lefts_idx[-1] = slice(0,-1)
    rights_idx[-1] = slice(1,None)
    lefts_idx, rights_idx, inner_idx = tuple(lefts_idx), tuple(rights_idx), tuple(inner_idx)

    # inner_idx here are where the interface fluxes are calculated with non-zero values.
    face_inner_idx = inner_idx
    u = np.zeros_like(Sol.rhoY)
    u[inner_idx] = 0.5 * (flux.rhoY[face_inner_idx][lefts_idx] + flux.rhoY[face_inner_idx][rights_idx]) / Sol.rhoY[inner_idx]

    shape = Sol.u.shape
    Diffs = States(shape,ud)
    Ampls = Characters(shape)
    Lefts = States(shape, ud)
    Rights = States(shape, ud)

    Diffs.u[...,:-1] = Sol.u[rights_idx] - Sol.u[lefts_idx]
    Diffs.v[...,:-1] = Sol.v[rights_idx] - Sol.v[lefts_idx]
    Diffs.w[...,:-1] = Sol.w[rights_idx] - Sol.w[lefts_idx]
    Diffs.X[...,:-1] = Sol.X[rights_idx] - Sol.X[lefts_idx]
    Diffs.Y[...,:-1] = 1.0 / Sol.Y[rights_idx] - 1.0 / Sol.Y[lefts_idx]

    Slopes = slopes(Sol, Diffs, ud, elem)

    Ampls.u[...] = 0.5 * Slopes.u * (1. - lmbda * u)
    Ampls.v[...] = 0.5 * Slopes.v * (1. - lmbda * u)
    Ampls.w[...] = 0.5 * Slopes.w * (1. - lmbda * u)
    Ampls.X[...] = 0.5 * Slopes.X * (1. - lmbda * u)
    Ampls.Y[...] = 0.5 * Slopes.Y * (1. - lmbda * u)
    
    Lefts.u[...] = Sol.u + order_two * Ampls.u
    Lefts.v[...] = Sol.v + order_two * Ampls.v
    Lefts.w[...] = Sol.w + order_two * Ampls.w
    Lefts.X[...] = Sol.X + order_two * Ampls.X
    Lefts.Y[...] = 1.0 / (1.0 / Sol.Y + order_two * Ampls.Y)

    Ampls.u[...] = -0.5 * Slopes.u * (1. + lmbda * u)
    Ampls.v[...] = -0.5 * Slopes.v * (1. + lmbda * u)
    Ampls.w[...] = -0.5 * Slopes.w * (1. + lmbda * u)
    Ampls.X[...] = -0.5 * Slopes.X * (1. + lmbda * u)
    Ampls.Y[...] = -0.5 * Slopes.Y * (1. + lmbda * u)

    Rights.u[...] = Sol.u + order_two * Ampls.u
    Rights.v[...] = Sol.v + order_two * Ampls.v
    Rights.w[...] = Sol.w + order_two * Ampls.w
    Rights.X[...] = Sol.X + order_two * Ampls.X
    Rights.Y[...] = 1.0 / (1.0 / Sol.Y + order_two * Ampls.Y)

    vel = [Sol.u,Sol.v,Sol.w]

    # Lefts.rhoY[lefts_idx] = Rights.rhoY[rights_idx] = 0.5 * (Sol.rhoY[lefts_idx] + Sol.rhoY[rights_idx]) \
    #     - order_two * 0.5 * lmbda * (Sol.u[rights_idx] * Sol.rhoY[rights_idx] - Sol.u[lefts_idx] * Sol.rhoY[lefts_idx])
    Lefts.rhoY[lefts_idx] = Rights.rhoY[rights_idx] = 0.5 * (Sol.rhoY[lefts_idx] + Sol.rhoY[rights_idx]) \
        - order_two * 0.5 * lmbda * (vel[split_step][rights_idx] * Sol.rhoY[rights_idx] - vel[split_step][lefts_idx] * Sol.rhoY[lefts_idx])

    Lefts.p0[lefts_idx] = Rights.p0[rights_idx] = Lefts.rhoY[lefts_idx]**gamm

    get_conservatives(Rights, ud, th)
    get_conservatives(Lefts, ud, th)

    # return Lefts, Rights, u, Diffs, Ampls, Slopes
    return Lefts, Rights

def slopes(Sol, Diffs, ud, elem):
    """
    Reconstruct the piecewise linear slopes in the cells from the piecewise constants in the cell and its neighbours.

    Parameters
    ----------
    Sol : :py:class:`management.variable.Vars`
        Solution data container
    Diffs : :py:class:`management.variable.States`
        Data container for the difference in the quantities between adjacent cells, (right - left).
    ud : :py:class:`inputs.user_data.UserDataInit`
        Data container for the initial conditions
    elem : :py:class:`discretization.kgrid.ElemSpaceDiscr`
        Data container for the cell grid.

    Returns
    -------
    :py:class:`management.variable.Characters`
        Reconstructed piecewise linear slopes in the cell.
    """ 
    limiter_type_velocity = ud.limiter_type_velocity
    limiter_type_scalar = ud.limiter_type_scalars

    ndim = elem.ndim
    lefts_idx, rights_idx = [slice(None,)] * ndim, [slice(None,)] * ndim
    lefts_idx[-1] = slice(0,-1)
    rights_idx[-1] = slice(1,None)
    lefts_idx, rights_idx = tuple(lefts_idx), tuple(rights_idx)

    # amplitudes of the state differences:
    # first lefts_idx removes the zero at the end
    # since differences always result in len-1
    # and the second indexing selects lefts and rights
    aul = Diffs.u[lefts_idx][lefts_idx]
    avl = Diffs.v[lefts_idx][lefts_idx]
    awl = Diffs.w[lefts_idx][lefts_idx]
    aXl = Diffs.X[lefts_idx][lefts_idx]
    aYl = Diffs.Y[lefts_idx][lefts_idx]

    aur = Diffs.u[lefts_idx][rights_idx]
    avr = Diffs.v[lefts_idx][rights_idx]
    awr = Diffs.w[lefts_idx][rights_idx]
    aXr = Diffs.X[lefts_idx][rights_idx]
    aYr = Diffs.Y[lefts_idx][rights_idx]

    Slopes = Characters(Diffs.u.shape)

    Slopes.u[...,1:-1] = limiters(limiter_type_velocity, aul, aur)
    Slopes.v[...,1:-1] = limiters(limiter_type_velocity, avl, avr)
    Slopes.w[...,1:-1] = limiters(limiter_type_velocity, awl, awr)
    Slopes.X[...,1:-1] = limiters(limiter_type_scalar, aXl, aXr)
    Slopes.Y[...,1:-1] = limiters(limiter_type_scalar, aYl, aYr)

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
    if limiter_type == LimiterType.NONE:
        return 0.5 * (al + ar)

def get_conservatives(U, ud, th):
    """
    Get advected (conservative) quantities at the left and right of the cell interfaces.

    Parameters
    ----------
    U : :py:class:`management.variable.States`
        `Lefts` and `Rights` corresponding to the values at the cell interfaces.
    ud : :py:class:`inputs.user_data.UserDataInit`
        Data container for the initial conditions
    th : :py:class:`physics.gas_dynamics.thermodynamic.ThermodynamicInit`
        Class container for the thermodynamical constants.
    """
    U.rho = U.rhoY / U.Y
    U.rhou = U.u * U.rho
    U.rhov = U.v * U.rho
    U.rhow = U.w * U.rho
    U.rhoY = U.Y * U.rho
    U.rhoX = U.X * U.rho
    
    sgn = np.sign(U.rhoY)
    p = sgn*np.abs(U.rhoY)**th.gamminv
    
    U.rhoe = rhoe(U.rho, U.u, U.v, U.w, p, ud, th)