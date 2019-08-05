import numpy as np
from scipy import signal

def recompute_advective_fluxes(flux, Sol):
    ################################################
    # 2D case for now - generalise in future
    ################################################
    
    rhoYu = Sol.rhoY * Sol.rhou.T / Sol.rho

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
    flux[1].rhoY[1:-1,1:-2] = signal.convolve2d(rhoYv[:-1,:],kernel_v, mode='valid').T / kernel_v.sum()