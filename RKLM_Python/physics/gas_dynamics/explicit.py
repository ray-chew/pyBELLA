import numpy as np
from inputs.boundary import set_explicit_boundary_data, check_flux_bcs
from physics.gas_dynamics.recovery import recovery
from physics.gas_dynamics.numerical_flux import hll_solver
from management.variable import Vars

from debug import find_nearest

class NO_OF_RK_STAGES(object):
    NO_OF_RK_STAGES = 3

class TimeIntegratorParams(NO_OF_RK_STAGES):
    def __init__(self):
        self.dt_frac = 0
        self.flux_frac = np.zeros((self.NO_OF_RK_STAGES,2))
        self.update_frac = np.zeros((self.NO_OF_RK_STAGES))
        self.multiD_updt = False


def advect(Sol, flux, dt, elem, odd, ud, th, mpv, writer = None):
    # double strang sweep
    time_step = 0.5 * dt
    # print(time_step)
    stage = 0
    if (odd):
        for split in range(elem.ndim):
            lmbda = time_step / elem.dx
            Sol.flip()
            # elem.flip()
            # flux[0], flux[1] = flux[1], flux[0]
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)
    else:
        for i_split in range(elem.ndim):
            split = elem.ndim - 1 - i_split
            lmbda = time_step / elem.dx
            print("Done - stage = 0, split = ", split)
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv, writer)
            Sol.flip()
            # elem.flip()
            # flux[0], flux[1] = flux[1], flux[0]

    stage = 1
    if (odd):
        for i_split in range(elem.ndim):
            split = elem.ndim - 1 - i_split
            lmbda = time_step / elem.dx
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)
            Sol.flip()
            # elem.flip()
            # flux[0], flux[1] = flux[1], flux[0]
    else:
        for split in range(elem.ndim):
            # i_split = elem.ndim - 1 - split
            lmbda = time_step / elem.dx
            Sol.flip()
            # elem.flip()
            print("Done - stage = 1, split = ", split)
            # flux[0], flux[1] = flux[1], flux[0]
            explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv, writer)
            
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

truefalse = True
counter = 0
def explicit_step_and_flux(Sol, flux, lmbda, elem, split_step, stage, ud, th, mpv, writer = None):
    set_explicit_boundary_data(Sol, elem, ud, th, mpv, dim=split_step)

    Lefts, Rights = recovery(Sol, flux, lmbda, ud, th, elem)

    global counter
    global truefalse
    if counter < 5 and writer != None:
        writer.populate('00' + str(counter),'Lefts_rhou',Lefts.rhou)
        counter += 1
    # if truefalse == True:
    #     writer.populate('000','Lefts_rhou',Lefts.rhou)
    #     truefalse = False

    # global truefalse
    # if truefalse == True:
    #     print(Lefts.u.shape)
    #     print(Lefts.u.flatten()[300])
    #     find_nearest(Lefts.u, 14.294938667740171)
    #     truefalse = False
        
    # global truefalse
    # if truefalse == True:
    #     print(stage, split_step)
    #     print(Lefts.u[0])
    #     # val, idx = find_nearest(Lefts.u,0.99999986475549874)

    #     truefalse = False

    # skipped check_flux_bcs for now; first debug other functions
    # check_flux_bcs(Lefts, Rights, elem, split_step, ud)

    # will need it for the test cases long waves and acoustic
    hll_solver(flux,Lefts,Rights,Sol, lmbda, ud, th)

    right_idx = (slice(None),slice(1,None))
    left_idx = (slice(None),slice(0,-1))
    # print(flux.rho[left_idx] - flux.rho[right_idx])
    Sol.rho += lmbda * (flux.rho[left_idx] - flux.rho[right_idx])
    Sol.rhou += lmbda * (flux.rhou[left_idx] - flux.rhou[right_idx])
    Sol.rhov += lmbda * (flux.rhov[left_idx] - flux.rhov[right_idx])
    Sol.rhow += lmbda * (flux.rhow[left_idx] - flux.rhow[right_idx])
    Sol.rhoe += lmbda * (flux.rhoe[left_idx] - flux.rhoe[right_idx])
    Sol.rhoY += lmbda * (flux.rhoY[left_idx] - flux.rhoY[right_idx])

    set_explicit_boundary_data(Sol, elem, ud, th, mpv, dim=split_step)



