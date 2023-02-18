# -*- coding: utf-8 -*-
import numpy as np
from scipy import signal
from management.debug import find_nearest

def recompute_advective_fluxes(flux, Sol, *args, **kwargs):
    """
    Recompute the advective fluxes at the cell interfaces, i.e. the faces. This function updates the `flux` container in-place.

    Parameters
    ----------
    flux : :py:class:`management.variable.States`
        Data container for the fluxes at the cell interfaces.
    Sol : :py:class:`management.variable.States`
        Data container for the Solution.

    Attention
    ---------
    This function is a mess and requires cleaning up.

    """
    ndim = Sol.rho.ndim
    inner_idx = tuple([slice(1,-1)] * ndim)

    if ndim == 2:
        kernel_u = np.array([[0.5,1.,0.5],[0.5,1.,0.5]])
        kernel_v = kernel_u.T
    elif ndim == 3:
        kernel_u = np.array([[[1,2,1],[2,4,2],[1,2,1]],[[1,2,1],[2,4,2],[1,2,1]]])
        kernel_v = np.swapaxes(kernel_u,1,0)
        kernel_w = np.swapaxes(kernel_u,2,0)

        rhoYw = Sol.rhoY * Sol.rhow / Sol.rho
        flux[2].rhoY[inner_idx] = signal.fftconvolve(rhoYw, kernel_w, mode='valid') / kernel_w.sum()
    else:
        assert(0, "Unsupported dimension in recompute_advective_flux")

    rhoYu = kwargs.get('u',Sol.rhoY * Sol.rhou / Sol.rho)

    flux[0].rhoY[inner_idx] = np.moveaxis(signal.fftconvolve(rhoYu, kernel_u, mode='valid') / kernel_u.sum(), 0, -1)

    rhoYv = kwargs.get('v',Sol.rhoY * Sol.rhov / Sol.rho)
    if ndim == 2:
        flux[1].rhoY[inner_idx] = signal.fftconvolve(rhoYv, kernel_v, mode='valid') / kernel_v.sum()
    elif ndim == 3:
        flux[1].rhoY[inner_idx] = np.moveaxis(signal.fftconvolve(rhoYv, kernel_v, mode='valid') / kernel_v.sum(), -1,0)
    # flux[1].rhoY[...,-1] = 0.

def hll_solver(flux, Lefts, Rights, Sol, lmbda, ud, th):
    """
    HLL solver for the Riemann problem. Chooses the advected quantities from `Lefts` or `Rights` based on the direction given by `flux`.

    Parameters
    ----------
    flux : :py:class:`management.variable.States`
        Data container for fluxes.
    Lefts : :py:class:`management.variable.States`
        Container for the quantities on the left of the cell interfaces.
    Rights : :py:class:`management.variable.States`
        Container for the quantities on the right of the cell interfaces.
    Sol : :py:class:`management.variable.Vars`
        Solution data container.
    lmbda : float
        :math:`\\frac{dt}{dx}`, where :math:`dx` is the grid-size in the direction of the substep.
    ud : :py:class:`inputs.user_data.UserDataInit`
        Class container for the initial condition.
    th : :py:class:`physics.gas_dynamics.thermodynamic.ThermodynamicInit`
        Class container for the thermodynamical constants.

    Returns
    -------
    :py:class:`management.variable.States`
        `flux` data container with the solution of the Riemann problem.
    """
    # flux: index 1 to end = Left[inner_idx]: index 0 to -1 = Right[inner_idx]: index 1 to end
    
    ndim = Sol.rho.ndim
    left_idx, right_idx, remove_cols_idx = [slice(None)] * ndim, [slice(None)] * ndim, [slice(None)] * ndim

    remove_cols_idx[-1] = slice(1,-1)
    left_idx[-1] = slice(0,-1)
    right_idx[-1] = slice(1,None)

    left_idx, right_idx, remove_cols_idx = tuple(left_idx), tuple(right_idx), tuple(remove_cols_idx)

    Lefts.primitives(th)
    Rights.primitives(th)
    
    upwind = 0.5 * (1.0 + np.sign(flux.rhoY))
    upl = upwind[right_idx]
    upr = (1.0 - upwind[left_idx]) 

    flux.rhou[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.u[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.u[right_idx])
    flux.rho[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * 1.0 + upr[right_idx] / Rights.Y[right_idx] * 1.0)

    flux.rhov[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.v[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.v[right_idx])
    flux.rhow[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.w[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.w[right_idx])
    flux.rhoX[remove_cols_idx] = flux.rhoY[remove_cols_idx] * (upl[left_idx] / Lefts.Y[left_idx] * Lefts.X[left_idx] + upr[right_idx] / Rights.Y[right_idx] * Rights.X[right_idx])

    return flux