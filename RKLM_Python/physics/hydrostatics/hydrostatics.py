import numpy as np

def hydrostatic_column(HydroState, HydroState_n, Y, Y_n, elem, node, th, ud):
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
    HydroState_n.rho0[:,igy] = rhoY0 / Y_n[:,igy]
    HydroState_n.rhoY0[:,igy] = rhoY0
    HydroState_n.Y0[:,igy] = Y_n[:,igy]
    HydroState_n.S0[:,igy] = 1.0 / Y_n[:,igy]
    HydroState_n.p0[:,igy] = p0
    HydroState_n.p20[:,igy] = pi0 / ud.Msq

    dys = np.array([-elem.dy] + [-elem.dy/2] + [elem.dy/2] + list(np.ones((icy-3)) * elem.dy))
    S_p = 1.0 / Y[:,:]
    S_m = np.zeros_like(S_p)
    S_m[:,1:3] = 1.0 / Y_n[:,igy].reshape(-1,1)
    # S_m[2] = Y_n[0,igy]
    S_m[:,0] = 1.0 / Y[:,1]
    S_m[:,3:] = 1.0 / Y[:,2:-1]

    # print('S_p + S_m = ', S_p + S_m)

    S_integral_p = dys * 0.5 * (S_p + S_m)
    S_integral_p[:,:igy] = np.cumsum(S_integral_p[:,:igy][:,::-1],axis=1)[:,::-1]
    S_integral_p[:,igy:] = np.cumsum(S_integral_p[:,igy:],axis=1)

    # print('S_int_p = ', S_integral_p)
    # print(S_integral_p)

    pi_hydro = pi0 - Gamma * g * S_integral_p
    p_hydro = pi_hydro**Gamma_inv
    rhoY_hydro = pi_hydro**gm1_inv

    HydroState.rho0[...] = rhoY_hydro * S_p
    HydroState.p0[...] = p_hydro
    HydroState.p20[...] = pi_hydro / ud.Msq
    HydroState.S0[...] = S_p
    HydroState.S10[...] = 0.0
    HydroState.Y0[...] = 1.0 / S_p
    HydroState.rhoY0[...] = rhoY_hydro

    # Sn_m = Sn_p
    Sn_p = 1.0 / Y[:,:]
    dys = np.ones((icy)) * elem.dy
    dys[:2] *= -1
    Sn_integral_p = dys * Sn_p
    Sn_integral_p[:,:2] = np.cumsum(Sn_integral_p[:,:2][:,::-1],axis=1)[:,::-1]
    Sn_integral_p[:,2:] = np.cumsum(Sn_integral_p[:,2:],axis=1)

    pi_hydro_n = pi0 - Gamma * g * Sn_integral_p
    rhoY_hydro_n = pi_hydro_n**gm1_inv
    
    HydroState_n.rhoY0[:,:igy] = rhoY_hydro_n[:,:igy]
    HydroState_n.Y0[:,:igy] = Y_n[:,:igy]
    HydroState_n.S0[:,:igy] = 1.0 / Y_n[:,:igy]
    HydroState_n.p0[:,:igy] = rhoY_hydro_n[:,:igy]**th.gamm
    HydroState_n.p20[:,:igy] = pi_hydro_n[:,:igy] / ud.Msq

    HydroState_n.rhoY0[:,igy:] = rhoY_hydro_n[:,igy:]
    HydroState_n.Y0[:,igy:] = Y_n[:,igy:]
    HydroState_n.S0[:,igy:] = 1.0 / Y_n[:,igy:]
    HydroState_n.p0[:,igy:] = rhoY_hydro_n[:,igy:]**th.gamm
    HydroState_n.p20[:,igy:] = pi_hydro_n[:,igy:] / ud.Msq

def hydrostatic_state(mpv, elem, node, th, ud):
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

np.set_printoptions(precision=18)

def hydrostatic_initial_pressure(Sol,mpv,elem,node,ud,th):
    Gammainv = th.Gammainv
    igy = node.igy
    dx = node.dx
    dy = node.dy

    beta = np.zeros((node.icx))
    bdpdx = np.zeros((node.icx))

    x_idx_m = slice(0,-1)
    x_idx_c = slice(1,None)
    y_idx = slice(igy,-igy)
    xn_idx = slice(1,-1)
    height = node.y[-igy-1]

    Pc = Sol.rhoY[x_idx_c,y_idx]
    # print(Pc.shape)
    Pm = Sol.rhoY[x_idx_m,y_idx]
    thc = Pc / Sol.rho[x_idx_c,y_idx]
    thm = Pm / Sol.rho[x_idx_m,y_idx]
    beta[xn_idx] = np.sum(0.5 * (Pm * thm + Pc * thc) * dy, axis=1)
    bdpdx[xn_idx] = np.sum(0.5 * (Pm * thm + Pc * thc) * (mpv.p2_cells[x_idx_c,y_idx] - mpv.p2_cells[x_idx_m,y_idx]) * dy, axis=1)

    beta *= Gammainv / height
    bdpdx *= Gammainv / height / dx
    print(bdpdx)
    