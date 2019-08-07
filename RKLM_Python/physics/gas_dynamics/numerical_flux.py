import numpy as np
from scipy import signal

def recompute_advective_fluxes(flux, Sol):
    ################################################
    # 2D case for now - generalise in future
    ################################################
    
    rhoYu = Sol.rhoY * Sol.rhou.T / Sol.rho

    # keep for debugging purpose... last checked, function is correct
    # will need to think of a better method to generalise from 2d

    # idx = 901
    # print("rhoYu_cm = ", rhoYu.flatten()[idx-52])
    # print("rhoYu_mm = ", rhoYu.flatten()[idx-52-1])
    # print("rhoYu_cc = ", rhoYu.flatten()[idx])
    # print("rhoYu_mc = ", rhoYu.flatten()[idx-1])
    # print("rhoYu_cp = ", rhoYu.flatten()[idx+52])
    # print("rhoYu_mp = ", rhoYu.flatten()[idx+52-1])

    kernel_u = np.array([[0.5, 0.5],[1., 1.],[0.5, 0.5]])
    flux[0].rhoY[1:-2,1:-1] = signal.convolve2d(rhoYu[:,:-1], kernel_u, mode='valid').T / kernel_u.sum()

    rhoYv = Sol.rhoY * Sol.rhov.T / Sol.rho
    kernel_v = np.array([[0.5,1.,0.5],[0.5,1.,0.5]])
    flux[1].rhoY[1:-1,1:-2] = signal.convolve2d(rhoYv[:-1,:], kernel_v, mode='valid').T / kernel_v.sum()


def hll_solver(flux, Lefts, Rights, Sol, lmbda, ud, th):

    # flux: index 1 to end = Left[inner_idx]: index 0 to -1 = Right[inner_idx]: index 1 to end
    inner_idx = (slice(1,-1),slice(1,-1))
    # flux_inner_idx = (slice(1,-2),slice(1,-1))
    flux_inner_idx = flux.flux_inner_idx
    left_idx = (slice(None),slice(0,-1))
    right_idx = (slice(None),slice(1,None))
    first_col = (slice(None),slice(0))

    Lefts.primitives(th)
    Rights.primitives(th)
    
    upwind = 0.5 * (1.0 + np.sign(flux.rhoY))
    
    # print(upwind.shape)
    # print(Lefts.u[0])
    upl = upwind / Lefts.Y[inner_idx]
    upr = (1.0 - upwind) / Rights.Y[inner_idx]
    print(flux.rhou.shape)
    
    flux.rhou[flux_inner_idx][:,1:] = flux.rhoY[right_idx] * (upl[right_idx] * Lefts.u[inner_idx][left_idx] + upr[right_idx] * Rights.u[inner_idx][right_idx])
    flux.rhou[first_col] = 0.0
    print(flux.rhou[1])
