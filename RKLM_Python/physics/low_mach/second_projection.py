from inputs.enum_bdry import BdryType
from inputs.boundary import set_explicit_boundary_data, set_ghostnodes_p2
from physics.low_mach.laplacian import stencil_9pt_2nd_try, stencil_9pt
from scipy import signal
import numpy as np
from itertools import product

from scipy.sparse.linalg import LinearOperator, bicgstab

import h5py

def euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, dt, alpha_diff):
    nc = node.sc
    rhs = np.zeros_like(mpv.p2_nodes)
    p2 = np.copy(mpv.p2_nodes[node.igx:-node.igx,node.igy:-node.igy])
    # p2 = np.copy(mpv.p2_nodes)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt)

    rhs[...], rhs_max = divergence_nodes(rhs,elem,node,Sol,mpv,ud)
    rhs /= dt

    if ud.is_compressible:
        rhs = rhs_from_p_old(rhs,node,mpv)

    # stencil_9pt_operator(elem,node,mpv,ud)
    # hf = h5py.File('/home/ray/git-projects/RKLM_Reference/RKLM_Reference/output_travelling_vortex_3d_48_with_initial_projection/low_Mach_gravity_comp/wplusx/wplusx_000.h5', 'r')
    # mpv.wplus[0][...] = (hf.get('Data-Set-2').value)
    # mpv.wplus[1][...] = (hf.get('Data-Set-2').value)
    # hf.close()

    # print(mpv.wcenter)

    lap2D = stencil_9pt_2nd_try(elem,node,mpv,ud)
    # lap2D = LinearOperator((53**2,53**2),matvec=lap2D)
    lap2D = LinearOperator((49**2,49**2),matvec=lap2D)
    # print(mpv.wplus[0])
    # rhs = rhs.T

    p2,_ = bicgstab(lap2D,rhs[node.igx:-node.igx,node.igy:-node.igy].reshape(-1,),x0=p2.reshape(-1,),maxiter=500)

    # p2,_ = bicgstab(lap2D,rhs.ravel(),x0=p2.ravel(),maxiter=200)
    # p2_full = p2.reshape(53,53)
    p2_full = np.zeros(nc).squeeze()
    p2_full[node.igx:-node.igx,node.igy:-node.igy] = p2.reshape(49,49)
    
    # hf = h5py.File('/home/ray/git-projects/RKLM_Reference/RKLM_Reference/output_travelling_vortex_3d_48_with_initial_projection/low_Mach_gravity_comp/pnew/pnew_000.h5', 'r')
    # p2_full = hf.get('Data-Set-2').value
    # hf.close()
    mpv.dp2_nodes[...] = np.copy(p2_full)
    correction_nodes(Sol,elem,node,mpv,p2_full,dt,1.0,ud)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    # mpv.dp2_nodes[...] = p2 - mpv.p2_nodes
    mpv.p2_nodes[...] = p2_full

    set_ghostnodes_p2(mpv.p2_nodes,node,ud)


def correction_nodes(Sol,elem,node,mpv,p,dt,crange,ud):
    ndim = node.ndim
    # coriolis = ud.coriolis_strength[0]
    time_offset = 3.0 - ud.acoustic_order
    vcorr = 1.

    igx,igy,igz = node.igs[0],node.igs[1],node.igs[2]
    dx = node.dx
    dy = node.dy
    oodx = 1.0/dx
    oody = 1.0/dy

    hplusx = mpv.wplus[0]
    hplusy = mpv.wplus[1]

    dSdy = (mpv.HydroState_n.S0[0,igy+1:-igy] - mpv.HydroState_n.S0[0,igy:-igy-1]) * oody
    inner_idx = (slice(igx,-igx),slice(igy,-igy))
    inner_eidx = (slice(igx,-igx-1),slice(igy,-igy-1))
    p = p[inner_idx]
    # print(p)

    indices = [idx for idx in product([slice(0,-1),slice(1,None)], repeat=ndim)]

    signs_x = [-1., -1., +1., +1.]
    signs_y = [-1., +1., -1., +1.]

    Dpx = 0
    Dpy = 0
    cnt = 0
    for index in indices:
        Dpx += signs_x[cnt] * p[index]
        Dpy += signs_y[cnt] * p[index]
        cnt += 1

    Dpx *= 0.5 * oodx
    Dpy *= 0.5 * oody

    thinv = Sol.rho[inner_idx] / Sol.rhoY[inner_idx]

    Sol.rhou[inner_idx] += -dt * thinv * hplusx[inner_eidx] * Dpx
    Sol.rhov[inner_idx] += -dt * thinv * hplusy[inner_eidx] * Dpy * vcorr
    Sol.rhoX[inner_idx] += -time_offset * dt * dSdy * Sol.rhov[inner_idx] * vcorr


def operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt):
    g = ud.gravity_strength[1]
    Msq = ud.Msq
    Gammainv = th.Gammainv

    ndim = node.ndim
    nonhydro = ud.nonhydrostasy
    dy = elem.dy
    icz = elem.icz

    time_offset = 3.0 - ud.acoustic_order

    coriolis = ud.coriolis_strength[0]

    ccenter = - (ud.compressibility * ud.Msq) * th.gm1inv / (dt**2) / time_offset
    cexp = 2.0 - th.gamm

    igs = elem.igs
    igx = igs[0]
    igy = igs[1]

    # # igz = igs[2]
    nindim = np.empty((ndim),dtype='object')
    innerdim = np.copy(nindim)
    eindim = np.empty((ndim),dtype='object')

    for dim in range(ndim):
        is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
        nindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
        innerdim[dim] = slice(igs[dim],-igs[dim])
        # eindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)
        eindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)

        if dim == 1:
            y_idx = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)
            right_idx = None if -igs[dim]+is_periodic == 0 else -igs[dim]+is_periodic
            y_idx1 = slice(igs[dim]-is_periodic+1, right_idx)

    strat = 2.0 * (mpv.HydroState_n.Y0[0,y_idx1] - mpv.HydroState_n.Y0[0,y_idx]) / (mpv.HydroState_n.Y0[0,y_idx1] + mpv.HydroState_n.Y0[0,y_idx]) / dy

    nindim = tuple(nindim)
    eindim = tuple(eindim)
    innerdim = tuple(innerdim)
    # nindim = (slice(None,),slice(None,))

    Y = Sol.rhoY[nindim] / Sol.rho[nindim]
    coeff = Gammainv * Sol.rhoY[nindim] * Y
    fsqsc = dt**2 * coriolis**2
    fimp = 1.0 / (1.0 + fsqsc)
    Nsqsc = time_offset * dt**2 * (g / Msq) * strat
    gimp = 1.0 / (nonhydro + Nsqsc)

    for dim in range(ndim):
        if dim == 1:
            mpv.wplus[dim][eindim] = coeff * gimp
        else:
            mpv.wplus[dim][eindim] = coeff * fimp

    # icx = elem.icx
    # icy = elem.icy

    # is_x_periodic = ud.bdry_type[0] == BdryType.PERIODIC
    # is_y_periodic = ud.bdry_type[1] == BdryType.PERIODIC

    # tmp_hplusx = np.zeros_like(mpv.wplus[0]).reshape(-1,)
    # tmp_hplusy = np.zeros_like(mpv.wplus[1]).reshape(-1,)
    # # coeff = coeff.ravel()

    # for j in range(igy - is_y_periodic, icy - igy + is_y_periodic):
    #     m = j * icx

    #     strat = 2.0 * (mpv.HydroState_n.Y0[0,j+1] - mpv.HydroState_n.Y0[0,j]) / (mpv.HydroState_n.Y0[0,j+1] + mpv.HydroState_n.Y0[0,j]) / dy

    #     for i in range(igx - is_x_periodic, icx - igx + is_x_periodic):
    #         n = m + i

    #         Y = Sol.rhoY.ravel()[n] / Sol.rho.ravel()[n]
    #         coeff = Gammainv * Sol.rhoY.ravel()[n] * Y
    #         fsqsc = dt**2 * coriolis**2
    #         fimp = 1.0 / (1.0 + fsqsc)
    #         Nsqsc = time_offset * dt**2 * (g / Msq) * strat
    #         gimp = 1.0 / (nonhydro + Nsqsc)

    #         tmp_hplusx[n] = coeff * fimp
    #         tmp_hplusy[n] = coeff * gimp

    # mpv.wplus[0][...] = tmp_hplusx.reshape(node.sc).squeeze()
    # mpv.wplus[1][...] = tmp_hplusy.reshape(node.sc).squeeze()

    kernel = (ndim - 1) * np.ones((2,2))

    for iz in range(icz):
        iz = slice(None,) if iz == 0 else iz

        mpv.wcenter[iz][innerdim] = ccenter * signal.convolve2d(Sol.rhoY[igx-1:-igx+1,igy-1:-igy+1]**cexp,kernel,mode='valid') / kernel.sum()

    scale_wall_node_values(mpv, node, ud)


def scale_wall_node_values(mpv, node, ud, factor=0.5):
    # Untested !!!!!!!!!!!
    ndim = node.ndim
    igs = node.igs

    wall_idx = np.empty((ndim), dtype=object)
    for dim in range(ndim):
        wall_idx[dim] = slice(igs[dim],-igs[dim])

    for dim in range(ndim):
        is_wall = ud.bdry_type[dim] == BdryType.WALL
        if is_wall:
            for direction in [-1,1]:
                wall_idx[dim] = (igs[dim] + 1) * direction
                mpv.wcenter[wall_idx] *= factor


def divergence_nodes(rhs,elem,node,Sol,mpv,ud):
    ndim = elem.ndim
    igs = elem.igs
    dxyz = elem.dxyz
    inner_idx = np.empty((ndim), dtype=object)

    for dim in range(ndim):
        is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC

        inner_idx[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
        
    indices = [idx for idx in product([slice(0,-1),slice(1,None)], repeat = ndim)]
    signs = [sgn for sgn in product([1,-1], repeat = ndim)]

    inner_idx = tuple(inner_idx)
    Sols = np.dstack((Sol.rhou[inner_idx], Sol.rhov[inner_idx], Sol.rhow[inner_idx]))

    oodxyz = 1./dxyz
    Y = Sol.rhoY[inner_idx] / Sol.rho[inner_idx]
    tmp_fxyz = 0.5**(ndim-1) * oodxyz[:ndim] * Sols[...,:ndim] * Y[...,None]
    
    count = 0
    for index in indices:
        rhs[inner_idx][index] += np.inner(signs[count], tmp_fxyz)
        count += 1

    rhs_max = np.amax(rhs[inner_idx]) if np.amax(rhs[inner_idx]) > 0 else 0
    return rhs, rhs_max


def rhs_from_p_old(rhs,node,mpv):
    igs = node.igs
    ndim = node.ndim

    assert ndim != 1, "Not implemented for 1D"
    
    inner_idx = np.empty((ndim), dtype=object)
    for dim in range(ndim):
        igs[dim] = slice(igs[dim],-igs[dim])

    rhs_hh = mpv.wcenter[inner_idx] * mpv.p2_nodes[inner_idx]
    rhs[inner_idx] += rhs_hh
    return rhs

