import numpy as np

def HydrostaticStates(mpv, elem, node, th, ud):
    Gamma = th.gm1 / th.gamm
    gamm = th.gamm
    gm1 = th.gm1
    Gamma_inv = 1.0 / Gamma
    gm1_inv = 1.0 / gm1

    icy = elem.icy
    igy = elem.igy

    rhoY0 = 1.0

    g = ud.gravity_strength[1]

    p0 = rhoY0**gamm
    pi0 = rhoY0**gm1
    mpv.HydroState_n.Y0[0,igy] = ud.stratification(0.0)
    mpv.HydroState_n.rhoY0[0,igy] = rhoY0
    mpv.HydroState_n.rho0[0,igy] = rhoY0 / ud.stratification(0.0)
    mpv.HydroState_n.S0[0,igy] = 1.0 / mpv.HydroState_n.Y0[0,igy]
    mpv.HydroState_n.p0[0,igy] = p0
    mpv.HydroState_n.p20[0,igy] = pi0 / ud.Msq

    y_p = elem.y[igy-1]
    S_p = 1.0 / ud.stratification(y_p)
    S_integral_p = -0.5 * elem.dy * 0.5 * (S_p + 1.0 / ud.stratification(0.0))

    yn_p = node.y[igy-1]
    Sn_p = 1.0 / ud.stratification(elem.y[igy-1])
    Sn_integral_p = -node.dy * Sn_p

    ###########################
    y_m = np.copy(y_p)
    y_p = y_m - elem.dy
    S_m = np.copy(S_p)
    S_p = 1.0 / ud.stratification(y_p)
    # S_integral_p -= np.arange(igy)[::-1] * elem.dy * 0.5 * (S_m + S_p)
    S_integral_p -= np.arange(igy) * elem.dy
    S_integral_p = S_integral_p[::-1]

    pi_hydro = pi0 - Gamma * g * S_integral_p
    p_hydro = pi_hydro**Gamma_inv
    rhoY_hydro = pi_hydro**gm1_inv

    mpv.HydroState.rhoY0[0,:igy] = rhoY_hydro
    mpv.HydroState.rho0[0,:igy] = rhoY_hydro * S_p
    mpv.HydroState.p0[0,:igy] = p_hydro
    mpv.HydroState.p20[0,:igy] = pi_hydro / ud.Msq
    mpv.HydroState.S0[0,:igy] = S_p
    mpv.HydroState.S10[0,:igy] = 0.0
    mpv.HydroState.Y0[0,:igy] = 1.0 / S_p

    yn_m = np.copy(yn_p)
    yn_p = yn_m - elem.dy
    # Sn_m = np.copy(Sn_p)
    Sn_p = 1.0 / ud.stratification(0.5*(yn_p + yn_m))
    Sn_integral_p -= np.arange(0,igy) * elem.dy
    Sn_integral_p = Sn_integral_p[::-1]

    pi_hydro_n = pi0 - Gamma * g * Sn_integral_p
    rhoY_hydro_n = pi_hydro_n**gm1_inv
    mpv.HydroState_n.rhoY0[0,:igy] = rhoY_hydro_n
    mpv.HydroState_n.Y0[0,:igy] = ud.stratification(0.5 * (y_p + y_m))
    mpv.HydroState_n.rho0[0,:igy] = rhoY_hydro_n / mpv.HydroState_n.Y0[0,:igy]
    mpv.HydroState_n.S0[0,:igy] = 1.0 / mpv.HydroState_n.Y0[0,:igy]
    mpv.HydroState_n.p0[0,:igy] = rhoY_hydro_n**th.gamm
    mpv.HydroState_n.p20[0,:igy] = pi_hydro_n / ud.Msq

    ###########################

    y_p = elem.y[igy]
    S_p = 1.0 / ud.stratification(y_p)
    S_integral_p = 0.5 * elem.dy * 0.5 * (S_p + 1.0 / ud.stratification(0.0))

    yn_p = node.y[igy+1]
    Sn_p = 1.0 / ud.stratification(y_p)
    Sn_integral_p = node.dy * Sn_p
    y_m = np.copy(y_p)
    y_p = y_m + elem.dy
    S_m = np.copy(S_p)
    S_p = 1.0 / ud.stratification(y_p)

    S_integral_p += np.arange(icy-igy)*elem.dy
    pi_hydro = pi0 - Gamma * g * S_integral_p
    p_hydro = pi_hydro**Gamma_inv
    rhoY_hydro = pi_hydro**gm1_inv

    mpv.HydroState.rho0[0,igy:] = rhoY_hydro * S_p
    mpv.HydroState.p0[0,igy:] = p_hydro
    mpv.HydroState.p20[0,igy:] = pi_hydro / ud.Msq
    mpv.HydroState.S0[0,igy:] = S_p
    mpv.HydroState.S10[0,igy:] = 0.0
    mpv.HydroState.Y0[0,igy:] = 1.0 / S_p
    mpv.HydroState.rhoY0[0,igy:] = rhoY_hydro

    yn_m = np.copy(yn_p)
    yn_p = yn_m + node.dy
    # Sn_m = Sn_p
    Sn_p = 1.0 / ud.stratification(0.5 * (yn_p + yn_m))
    Sn_integral_p += np.arange(icy-igy)*elem.dy

    pi_hydro_n = pi0 - Gamma * g * Sn_integral_p
    rhoY_hydro_n = pi_hydro_n**gm1_inv

    mpv.HydroState_n.rhoY0[0,igy+1:] = rhoY_hydro_n
    mpv.HydroState_n.Y0[0,igy+1:] = ud.stratification(0.5 * (y_p + y_m))
    mpv.HydroState_n.rho0[0,igy+1:] = rhoY_hydro_n / mpv.HydroState_n.Y0[0,igy+1:]
    mpv.HydroState_n.S0[0,igy+1:] = 1.0 / mpv.HydroState_n.Y0[0,igy+1:]
    mpv.HydroState_n.p0[0,igy+1:] = rhoY_hydro_n**th.gamm
    mpv.HydroState_n.p20[0,igy+1:] = pi_hydro_n / ud.Msq

