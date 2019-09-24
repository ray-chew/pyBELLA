import numpy as np

from inputs.boundary import set_ghostcells_p2

def nonhydrostasy(ud,t):
    if ud.is_nonhydrostatic == 0:
        return 0.0
    elif ud.is_nonhydrostatic == 1:
        return 1.0
    elif ud.is_nonhydrostatic == -1:
        a = 12.5
        b = 1.0 / 24.0
        c = np.min(1.0, np.max(0.0, b * (t-a)))
        return c
    else:
        assert 0

def compressibility(ud,t,step):
    if ud.is_compressible == 0:
        return 0.0
    elif ud.is_compressible == 1:
        return 1.0
    elif ud.is_compressible == -1:
        # dtloc = 12.0
        # a = 0.0 * dtloc
        # b = 1.0 / dtloc / 20.0
        # tau = np.min(1.0, np.max(0.0 * b * (t-a)))
        # c = 0.5 * (1.0 - np.cos(np.pi * tau))
        # return c
        current_transition_step = step - ud.no_of_initial
        print(current_transition_step)
        print(np.linspace(0.0,1.0,ud.no_of_transition+2)[1:-1][current_transition_step])
        return np.linspace(0.0,1.0,ud.no_of_transition+2)[1:-1][current_transition_step]
    else:
        assert 0

def is_compressible(ud,step):
    if ud.continuous_blending == True:
        if step < ud.no_of_initial:
            return 0
        elif step < (ud.no_of_initial + ud.no_of_transition):
            return -1
        else:
            return 1
    else:
        return ud.is_compressible


def rhoe(rho,u,v,w,p,ud,th):
    Msq = ud.compressibility * ud.Msq
    gm1inv = th.gm1inv

    return p * gm1inv + 0.5 * Msq * rho * (u**2 + v**2 + w**2)

def synchronise_variables(mpv, Sol, elem, node, ud, th):
    scale_factor = 1.0 / ud.Msq

    if (ud.is_compressible):
        p2bg = mpv.HydroState.p20[0,:].reshape(1,-1)
        p2bg = np.repeat(p2bg, elem.icx,axis=0)
        mpv.p2_cells = scale_factor * Sol.rhoY**th.gm1 - p2bg

    set_ghostcells_p2(mpv.p2_cells, elem, ud)