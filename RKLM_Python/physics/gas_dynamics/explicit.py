import numpy as np
from inputs.boundary import set_explicit_boundary_data, check_flux_bcs
from .recovery import recovery
from .numerical_flux import hll_solver
from management.variable import Vars

def advect(Sol, flux, dt, elem, odd, ud, th, mpv, writer = None):
    # double strang sweep
    time_step = 0.5 * dt
    # print(time_step)
    ndim = elem.ndim
    stage = 0

    if (odd):
        for split in range(ndim):
            lmbda = time_step / elem.dxyz[split]
            Sol.flip_forward()
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)
    else:
        for i_split in range(ndim):
            split = elem.ndim - 1 - i_split
            lmbda = time_step / elem.dxyz[split]
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv, writer)
            Sol.flip_backward()

    stage = 1
    if (odd):
        for i_split in range(ndim):
            split = elem.ndim - 1 - i_split
            lmbda = time_step / elem.dxyz[split]
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)
            Sol.flip_backward()
    else:
        for split in range(ndim):
            lmbda = time_step / elem.dxyz[split]
            Sol.flip_forward()
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv, writer)
            
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

def explicit_step_and_flux(Sol, flux, lmbda, elem, split_step, stage, ud, th, mpv, writer = None):
    set_explicit_boundary_data(Sol, elem, ud, th, mpv, step=split_step)

    Lefts, Rights = recovery(Sol, flux, lmbda, ud, th, elem)

    # skipped check_flux_bcs for now; first debug other functions
    # check_flux_bcs(Lefts, Rights, elem, split_step, ud)

    # will need it for the test cases long waves and acoustic
    hll_solver(flux,Lefts,Rights,Sol, lmbda, ud, th)

    ndim = elem.ndim
    left_idx, right_idx = [slice(None)] * ndim, [slice(None)] * ndim
    right_idx[-1] = slice(1,None)
    left_idx[-1] = slice(0,-1)
    left_idx, right_idx = tuple(left_idx), tuple(right_idx)

    Sol.rho += lmbda * (flux.rho[left_idx] - flux.rho[right_idx])
    Sol.rhou += lmbda * (flux.rhou[left_idx] - flux.rhou[right_idx])
    Sol.rhov += lmbda * (flux.rhov[left_idx] - flux.rhov[right_idx])
    Sol.rhow += lmbda * (flux.rhow[left_idx] - flux.rhow[right_idx])
    Sol.rhoe += lmbda * (flux.rhoe[left_idx] - flux.rhoe[right_idx])
    Sol.rhoX += lmbda * (flux.rhoX[left_idx] - flux.rhoX[right_idx])
    Sol.rhoY += lmbda * (flux.rhoY[left_idx] - flux.rhoY[right_idx])
    
    set_explicit_boundary_data(Sol, elem, ud, th, mpv, step=split_step)