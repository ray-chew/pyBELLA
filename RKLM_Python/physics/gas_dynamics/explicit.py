import numpy as np
from inputs.boundary import set_explicit_boundary_data
from physics.gas_dynamics.recovery import recovery
from physics.gas_dynamics.numerical_flux import hll_solver
from management.variable import Vars

class NO_OF_RK_STAGES(object):
    NO_OF_RK_STAGES = 3

class TimeIntegratorParams(NO_OF_RK_STAGES):
    def __init__(self):
        self.dt_frac = 0
        self.flux_frac = np.zeros((self.NO_OF_RK_STAGES,2))
        self.update_frac = np.zeros((self.NO_OF_RK_STAGES))
        self.multiD_updt = False


def advect(Sol, flux, dt, elem, odd, ud, th, mpv):
    # double strang sweep
    time_step = 0.5 * dt

    stage = 0
    if (odd):
        for split in range(elem.ndim):
            lmbda = time_step / elem.dx
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)
    else:
        for i_split in range(elem.ndim):
            split = elem.ndim - 1 - i_split
            Sol.flip()
            flux[split].trim_zeros()
            flux[split].flip()
            flux[split].get_flux_inner_idx(i_split)
            lmbda = time_step / elem.dx
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)

    stage = 1
    if (odd):
        for i_split in range(elem.ndim):
            split = elem.ndim - 1 - i_split
            lmbda = time_step / elem.dx
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)
    else:
        for split in range(elem.ndim):
            lmbda = time_step / elem.dx
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)


def explicit_step_and_flux(Sol, flux, lmbda, elem, split_step, stage, ud, th, mpv):
    set_explicit_boundary_data(Sol, elem, ud, th, mpv, dim = split_step)

    Lefts, Rights = recovery(Sol, flux, lmbda, ud, th, elem)
    # skipped check_flux_bcs for now; first debug other functions
    # will need it for the test cases long waves and acoustic
    hll_solver(flux,Lefts,Rights,Sol, lmbda, ud, th)

