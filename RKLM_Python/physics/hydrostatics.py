from inputs.boundary import set_ghostcells_p2, set_ghostnodes_p2
import numpy as np
import numba

def hydrostatic_column(HydroState, HydroState_n, Y, Y_n, elem, node, th, ud):
    Gamma = th.gm1 / th.gamm
    gamm = th.gamm
    gm1 = th.gm1
    Gamma_inv = 1.0 / Gamma
    gm1_inv = 1.0 / gm1

    icy = elem.icy
    igy = elem.igy

    xc_idx = slice(0,-1)
    yc_idx = slice(0,-1)

    c_idx = (xc_idx,yc_idx)

    rhoY0 = 1.0

    g = ud.gravity_strength[1]

    p0 = rhoY0**gamm
    pi0 = rhoY0**gm1
    HydroState_n.rho0[xc_idx,igy] = rhoY0 / Y_n[:,igy]
    HydroState_n.rhoY0[xc_idx,igy] = rhoY0
    HydroState_n.Y0[xc_idx,igy] = Y_n[:,igy]
    HydroState_n.S0[xc_idx,igy] = 1.0 / Y_n[:,igy]
    HydroState_n.p0[xc_idx,igy] = p0
    HydroState_n.p20[xc_idx,igy] = pi0 / ud.Msq

    dys = np.array([-elem.dy] + [-elem.dy/2] + [elem.dy/2] + list(np.ones((icy-3)) * elem.dy))
    S_p = 1.0 / Y[:,:]
    S_m = np.zeros_like(S_p)
    S_m[:,igy-1:igy+1] = 1.0 / Y_n[:,igy].reshape(-1,1)
    S_m[:,0] = 1.0 / Y[:,igy-1]
    S_m[:,igy+1:] = 1.0 / Y[:,igy:-1]

    S_integral_p = dys * 0.5 * (S_p + S_m)
    S_integral_p[:,:igy] = np.cumsum(S_integral_p[:,:igy][:,::-1],axis=1)[:,::-1]
    S_integral_p[:,igy:] = np.cumsum(S_integral_p[:,igy:],axis=1)

    pi_hydro = pi0 - Gamma * g * S_integral_p
    p_hydro = pi_hydro**Gamma_inv
    rhoY_hydro = pi_hydro**gm1_inv

    HydroState.rho0[c_idx] = rhoY_hydro * S_p
    HydroState.p0[c_idx] = p_hydro
    HydroState.p20[c_idx] = pi_hydro / ud.Msq
    HydroState.S0[c_idx] = S_p
    HydroState.S10[c_idx] = 0.0
    HydroState.Y0[c_idx] = 1.0 / S_p
    HydroState.rhoY0[c_idx] = rhoY_hydro

    Sn_p = 1.0 / Y[:,:]
    dys = np.ones((icy)) * elem.dy
    dys[:igy] *= -1
    Sn_integral_p = dys * Sn_p
    Sn_integral_p[:,:igy] = np.cumsum(Sn_integral_p[:,:igy][:,::-1],axis=1)[:,::-1]
    Sn_integral_p[:,igy:] = np.cumsum(Sn_integral_p[:,igy:],axis=1)

    pi_hydro_n = pi0 - Gamma * g * Sn_integral_p
    rhoY_hydro_n = pi_hydro_n**gm1_inv
    
    HydroState_n.rhoY0[xc_idx,:igy] = rhoY_hydro_n[:,:igy]
    HydroState_n.Y0[xc_idx,:igy] = Y_n[0,:igy]
    HydroState_n.S0[xc_idx,:igy] = 1.0 / Y_n[:,:igy]
    HydroState_n.p0[xc_idx,:igy] = rhoY_hydro_n[:,:igy]**th.gamm
    HydroState_n.p20[xc_idx,:igy] = pi_hydro_n[:,:igy] / ud.Msq

    HydroState_n.rhoY0[xc_idx,igy+1:] = rhoY_hydro_n[:,igy:]
    HydroState_n.Y0[xc_idx,igy+1:] = Y_n[0,igy:]
    HydroState_n.S0[xc_idx,igy+1:] = 1.0 / Y_n[:,igy:]
    HydroState_n.p0[xc_idx,igy+1:] = rhoY_hydro_n[:,igy:]**th.gamm
    HydroState_n.p20[xc_idx,igy+1:] = pi_hydro_n[:,igy:] / ud.Msq


def hydrostatic_state(mpv, elem, node, th, ud):
    g = ud.gravity_strength[1]
    Gamma = th.Gamma
    Hex = 1.0 / (th.Gamma * g)
    dy = elem.dy

    pi_np = np.exp(-(node.y + 0.5 * dy) / Hex)
    pi_nm = np.exp(-(node.y - 0.5 * dy) / Hex)
    pi_n = np.exp(-(node.y) / Hex)

    Y_n = - Gamma * g * dy / (pi_np - pi_nm)
    P_n = pi_n**th.gm1inv
    p_n = pi_n**th.Gammainv
    rho_n = P_n / Y_n

    mpv.HydroState_n.p20[...] = pi_n / ud.Msq
    mpv.HydroState_n.p0[...] = p_n
    mpv.HydroState_n.rho0[...] = rho_n
    mpv.HydroState_n.rhoY0[...] = P_n
    mpv.HydroState_n.Y0[...] = Y_n
    mpv.HydroState_n.S0[...] = 1.0 / Y_n

    pi_cp = np.exp(-(elem.y + 0.5 * dy) / Hex)
    pi_cm = np.exp(-(elem.y - 0.5 * dy) / Hex)
    pi_c  = np.exp(-(elem.y) / Hex)

    Y_c = - Gamma * g * dy / (pi_cp - pi_cm)
    P_c = pi_c**th.gm1inv
    p_c = pi_c**th.Gammainv
    rho_c = P_c / Y_c

    mpv.HydroState.p20[...] = pi_c / ud.Msq
    mpv.HydroState.p0[...] = p_c
    mpv.HydroState.rho0[...] = rho_c
    mpv.HydroState.rhoY0[...] = P_c
    mpv.HydroState.Y0[...] = Y_c
    mpv.HydroState.S0[...] = 1.0 / Y_c

    

def hydrostatic_initial_pressure(Sol,mpv,elem,node,ud,th):
    Gammainv = th.Gammainv
    igy = node.igy
    igx = node.igx
    icx= elem.icx
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
    Pm = Sol.rhoY[x_idx_m,y_idx]
    thc = Pc / Sol.rho[x_idx_c,y_idx]
    thm = Pm / Sol.rho[x_idx_m,y_idx]
    beta[xn_idx] = np.sum(0.5 * (Pm * thm + Pc * thc) * dy, axis=1)
    bdpdx[xn_idx] = np.sum(0.5 * (Pm * thm + Pc * thc) * (mpv.p2_cells[x_idx_c,y_idx] - mpv.p2_cells[x_idx_m,y_idx]) * dy, axis=1)

    beta *= Gammainv / height
    bdpdx *= Gammainv / height / dx

    coeff = np.zeros((elem.icx))
    pibot = np.zeros((elem.icx))
    coeff[igx+1:-igx+1] = np.cumsum(coeff[igx:-igx] + dx / beta[igx+1:-igx])
    pibot[igx+1:-igx+1] = np.cumsum(pibot[igx:-igx] - dx * bdpdx[igx+1:-igx] / beta[igx+1:-igx])
    
    dotPU = pibot[icx-igx] / coeff[icx-igx]
    pibot[igx:-igx] -= dotPU * coeff[igx:-igx]

    x_idx = slice(igx,-igx+1)
    y_idx = slice(igy,-igy+1)

    mpv.p2_cells[x_idx,y_idx] += pibot[x_idx].reshape(-1,1) - 1.0 * mpv.HydroState.p20[y_idx].reshape(1,-1)
    set_ghostcells_p2(mpv.p2_cells, elem, ud)

    icxn = node.icx
    icyn = node.icy
    x_idx = slice(1,icxn-1)
    y_idx = slice(igy,-igy+1)
    height = node.y[-igy]

    Pc = Sol.rhoY[1:,y_idx]
    thc = Pc / Sol.rho[1:,y_idx]

    beta = np.zeros((elem.icx,))
    bdpdx = np.zeros((elem.icx))
    
    beta[1:] = np.sum(Pc*thc*dy, axis=1)
    beta *= Gammainv / height

    bdpdx[1:] = np.sum(Pc * thc * (mpv.p2_nodes[1:-1,igy:-igy] - mpv.p2_nodes[:-2,igy:-igy]) * dy, axis=1)
    bdpdx *= Gammainv / height / dx

    coeff = np.zeros((node.icx))
    pibot = np.zeros((node.icx))

    coeff[igx+1:-igx+1] = np.cumsum(coeff[igx:-igx] + dx / beta[igx+1:])
    pibot[igx+1:-igx+1] = np.cumsum(pibot[igx:-igx] - dx * bdpdx[igx+1:] / beta[igx+1:])    

    dotPU = pibot[icx-igx] / coeff[icx-igx]

    pibot[igx:-igx] -= dotPU * coeff[igx:-igx]

    x_idx = slice(igx,-igx+1)
    y_idx = slice(igy,-igy+1)
    mpv.p2_nodes[x_idx,y_idx] += pibot[x_idx].reshape(-1,1) - 1.0 * mpv.HydroState_n.p20[y_idx].reshape(1,-1)

    mpv.dp2_nodes[:,:] = mpv.p2_nodes

    # guess initial node value (at left-most node)
    mpv.p2_nodes[igx,igy:-igy] = mpv.dp2_nodes[igx,igy:-igy]

    mpv.p2_nodes[:,:] = loop_over_array(igx,igy,icxn,icyn,mpv.p2_nodes, mpv.dp2_nodes)

    assert ((node.icx+1)%2) == 1
    delp2 = 0.5 * (mpv.p2_nodes[-igx-1,igy:-igy] - mpv.p2_nodes[igx,igy:-igy])
    delp2 = delp2.reshape(1,-1)
    sgn = np.ones_like(mpv.p2_nodes[:,0][igy:-igy]).reshape(-1,1)
    
    sgn[1::2] *= -1

    mpv.p2_nodes[igx:-igx,igy:-igy] += sgn * delp2
    set_ghostnodes_p2(mpv.p2_nodes, node, ud)

    mpv.dp2_nodes[:,:] = 0.0

    inner_domain = (slice(igx,-igx), slice(igy,-igy))
    pi = ud.Msq * (mpv.p2_cells[inner_domain] + 1.0 * mpv.HydroState.p20[igy:-igy])
    Y = Sol.rhoY[inner_domain] / Sol.rho[inner_domain]
    rhoold = np.copy(Sol.rho[inner_domain])
    Sol.rhoY[inner_domain] = pi**th.gm1inv
    Sol.rho[inner_domain] = Sol.rhoY[inner_domain] / Y
    Sol.rhou[inner_domain] *= Sol.rho[inner_domain] / rhoold
    Sol.rhov[inner_domain] *= Sol.rho[inner_domain] / rhoold
    Sol.rhow[inner_domain] *= Sol.rho[inner_domain] / rhoold
    Sol.rhoX[inner_domain] *= Sol.rho[inner_domain] / rhoold

# need details:
# populate the rest of the nodes recursively based on the left-most node.
# recursive: use numba.
@numba.jit(nopython=True)
def loop_over_array(igx,igy,icxn,icyn,p,dp):
    for j in range(igy,icyn-igy):
        for i in range(igx+1,icxn-igx):
            p[i,j] = 2.0 * dp[i-1,j] - p[i-1,j]
    return p
        

    