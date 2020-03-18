import numpy as np

from inputs.boundary import set_ghostcells_p2

def nonhydrostasy(ud,t,step):
    if step >= 0:
        if ud.is_nonhydrostatic == 0:
            return 0.0
        elif ud.is_nonhydrostatic == 1:
            return 1.0
        elif ud.is_nonhydrostatic == -1:
            current_transition_step = step - ud.no_of_hy_initial
            # print("current_transition_step =", step - ud.no_of_hy_initial)
            # print(np.linspace(0.0,1.0,ud.no_of_hy_transition+2)[1:-1][current_transition_step])
            return np.linspace(0.0,1.0,ud.no_of_hy_transition+2)[1:-1][current_transition_step]
    else:
        return float(ud.is_nonhydrostatic)

def compressibility(ud,t,step):
    if step >= 0:
        if ud.is_compressible == 0:
            if step < ud.no_of_pi_initial:
                return 0.0
            if ud.continuous_blending == False:
                return 0.0
            else:
                current_transition_step = step - ud.no_of_pi_initial
                return np.linspace(0.0,1.0,ud.no_of_pi_transition+2)[1:-1][current_transition_step]
        elif ud.is_compressible == 1:
            return 1.0
    # elif ud.is_compressible == 0:
    else:
        return ud.compressibility

def is_compressible(ud,step):
    if step >= 0:
        if ud.continuous_blending == True:
            if step < ud.no_of_pi_initial:
                return 0
            elif step < (ud.no_of_pi_initial + ud.no_of_pi_transition):
                return 0
            else:
                return 1
    else:
        return ud.is_compressible

def is_nonhydrostatic(ud,step):
    # print("is_nonhydrostatic", ud.is_nonhydrostatic)
    # print("no_of_nhy_initial:", ud.no_of_hy_initial)
    # print("no_of_nhy_transition:", ud.no_of_hy_transition)
    if step >= 0:
        if ud.continuous_blending == True:
            if step < ud.no_of_hy_initial:
                return 0
            elif step < (ud.no_of_hy_initial + ud.no_of_hy_transition):
                return -1
            else:
                return 1
    else:
        return ud.is_nonhydrostatic


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