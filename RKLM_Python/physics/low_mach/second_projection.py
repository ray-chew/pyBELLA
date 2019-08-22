from inputs.enum_bdry import BdryType
from inputs.boundary import set_explicit_boundary_data, set_ghostnodes_p2
from physics.low_mach.laplacian import stencil_9pt_3rd_try, stencil_9pt_2nd_try, stencil_9pt
from scipy import signal
import numpy as np
from itertools import product

from scipy.sparse.linalg import LinearOperator, bicgstab, gmres
from scipy.sparse import eye
import matplotlib.pyplot as plt

from debug import find_nearest
import h5py

def euler_forward_non_advective(Sol, mpv, elem, node, dt, ud, th):
    nonhydro = ud.nonhydrostasy
    wp = 1. # with pressure

    g = ud.gravity_strength[1]
    Msq = ud.Msq
    Ginv = th.Gammainv
    coriolis = ud.coriolis_strength[0]
    u0 = ud.wind_speed

    p2n = np.copy(mpv.p2_nodes)
    dp2n = np.zeros_like(p2n.T)
    dx = node.dy
    dy = node.dx

    mpv.rhs = np.array(mpv.rhs)
    mpv.rhs *= 0.0
    mpv.rhs[...], _ = divergence_nodes(mpv.rhs,elem,node,Sol,ud)
    div = mpv.rhs.T

    ######## Test this! #########
    scale_wall_node_values(mpv.rhs, node, ud, 2.0)
    #############################

    ## 2D-case ###
    inner_idx = (slice(1,-1),slice(1,-1))
    p2n = p2n[inner_idx]

    dpdx_kernel = np.array([[-1.,1.],[-1.,1.]])
    dpdy_kernel = np.array([[-1.,-1.],[1.,1.]])

    dpdx = -wp * 0.5 * signal.convolve2d(p2n, dpdx_kernel, mode='valid') / dx
    dpdy = -wp * 0.5 * signal.convolve2d(p2n, dpdy_kernel, mode='valid') / dy

    dpdx = dpdx.T
    dpdy = dpdy.T
    # div = div.T
    # if Sol.rhou.shape[0] == Sol.rhou.shape[1]:
    #     Sol.rhou[...] = Sol.rhou.T
    #     Sol.rhov[...] = Sol.rhov.T
    #     Sol.rhow[...] = Sol.rhow.T
    Sol.flip()
    # Sol.rhou = Sol.rhou.T
    # Sol.rhov = Sol.rhov.T
    # Sol.rhow = Sol.rhow.T
    # Sol.rho = Sol.rho.T
    # Sol.rhoY = Sol.rhoY.T

    # idx = 587
    # print('Sol.rhou(before) = ', Sol.rhou.T.flatten()[idx])
    # print('Sol.rhov(before) = ', Sol.rhov.T.flatten()[idx])
    # find_nearest(Sol.rhou,0.50010834727366982)

    rhoYovG = Ginv * Sol.rhoY[inner_idx]
    dchi = 0. ###### INCOMPLETE!
    dbuoy = -1. * Sol.rhoY[inner_idx] * dchi
    drhou = Sol.rhou[inner_idx] - u0 * Sol.rho[inner_idx]

    rhoY = Sol.rhoY**(th.gamm - 2.0)    
    dpidP_kernel = np.array([[1.,1.],[1.,1.]])
    dpidP = (th.gm1 / ud.Msq) * signal.convolve2d(rhoY, dpidP_kernel, mode='valid') / dpidP_kernel.sum()
    
    Sol.rhou[inner_idx] = Sol.rhou[inner_idx] + dt * ( -1. * rhoYovG * dpdx + coriolis * Sol.rhow[inner_idx])
    Sol.rhov[inner_idx] = Sol.rhov[inner_idx] + dt * ( -1. * rhoYovG * dpdy + (g/Msq) * dbuoy) * nonhydro
    Sol.rhow[inner_idx] = Sol.rhow[inner_idx] - dt * coriolis * drhou

    # idx = 585
    # n_rows = int(np.floor(idx / 52))
    # n_cols = idx - (n_rows * 52)
    # n_rows -= 1
    # n_cols -= 1
    # idx = n_rows * 50 + n_cols

    # idx_n = 598
    # n_rows = int(np.floor(idx_n / 53))
    # n_cols = idx_n - (n_rows * 53) 
    # n_rows -= 1
    # n_cols -= 1
    # idx_n = n_rows * 51 + n_cols

    # print(idx_n)

    # print('dt = %.8f, dx = %.8f, dy = %.8f' %(dt,dx,dy))
    # print('p2n[idx_n] = ', p2n.flatten()[idx_n])
    # print('dpdx = ', dpdx.flatten()[idx])
    # print('dpdy = ', dpdy.flatten()[idx])
    # print('rhoYovG = ', rhoYovG.flatten()[idx])
    # print('Sol.rhou(after) = ', Sol.rhou[inner_idx].flatten()[idx])
    # print('Sol.rhov(after) = ', Sol.rhov[inner_idx].flatten()[idx])

    dp2n[inner_idx] -= dt * dpidP * div[inner_idx]
    # find_nearest(dp2n,-5.5631383332701323e-07)

    # idx_dp2n = 758
    # print('dp2n = ', dp2n.flatten()[idx_dp2n])
    if (ud.is_compressible):
        weight = ud.acoustic_order - 1.0
        mpv.p2_nodes += weight * dp2n.T

    # if Sol.rhou.shape[0] == Sol.rhou.shape[1]:
    #     Sol.rhou[...] = Sol.rhou.T
    #     Sol.rhov[...] = Sol.rhov.T
    #     Sol.rhow[...] = Sol.rhow.T
    Sol.flip()

    set_ghostnodes_p2(mpv.p2_nodes,node, ud)
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)
    
    # Sol.rhou = Sol.rhou.T
    # Sol.rhov = Sol.rhov.T
    # Sol.rhow = Sol.rhow.T
    # Sol.rho = Sol.rho.T
    # Sol.rhoY = Sol.rhoY.T

def euler_backward_non_advective_expl_part(Sol, mpv, elem, dt, ud, th):
    nonhydro = ud.nonhydrostasy

    # y-axis = gravity-direction
    g = ud.gravity_strength[1]
    Msq = ud.Msq
    dy = elem.dy

    time_offset = 3.0 - ud.acoustic_order
    coriolis = ud.coriolis_strength[0]
    u0 = ud.wind_speed
    fsqsc = dt**2 * coriolis**2
    ooopfsqsc = 1.0 / (1.0 + fsqsc)

    first_nodes_row_right_idx = (slice(0,1), slice(1,None))
    first_nodes_row_left_idx = (slice(0,1),slice(0,-1))

    strat = 2.0 * (mpv.HydroState_n.Y0[first_nodes_row_right_idx] - mpv.HydroState_n.Y0[first_nodes_row_left_idx]) / (mpv.HydroState_n.Y0[first_nodes_row_right_idx] + mpv.HydroState_n.Y0[first_nodes_row_left_idx]) / dy
    # print(mpv.HydroState_n.Y0[0])
    # works for 2d grid - check for non-square
    strat = strat.reshape(1,-1)
    strat = np.repeat(strat, elem.icx, axis=0)
    # print(strat)

    Nsqsc = time_offset * dt**2 * (g / Msq) * strat
    dbouy = 0.0 ###### Incomplete!!! ######
    rhov = (nonhydro * Sol.rhov + dt * (g/Msq) * dbouy) / (nonhydro + Nsqsc)

    drhou = Sol.rhou - u0 * Sol.rho
    Sol.rhou[...] = u0 * Sol.rho + ooopfsqsc * (drhou + dt * coriolis * Sol.rhow)
    Sol.rhov[...] = rhov
    Sol.rhow[...] = ooopfsqsc * (Sol.rhow - dt * coriolis * drhou)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)



def euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, dt, alpha_diff, writer = None):
    nc = node.sc
    rhs = np.zeros_like(mpv.p2_nodes)

    p2 = np.copy(mpv.p2_nodes[node.igx:-node.igx,node.igy:-node.igy])
    
    if writer != None:
        writer.populate('after_ebnaimp','p2_initial',mpv.p2_nodes)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt)
    
    base_filename = '/home/ray/git-projects/RKLM_Reference/RKLM_Reference/output_acoustic_wave_high/low_Mach_gravity_comp/'
    # p2 = h5py.File(base_filename + 'p2_initial/p2_initial_004.h5', 'r')['Data-Set-2'][:]
    # p2 = p2[2:-2,2:-2]
    # mpv.wcenter[...] = h5py.File(base_filename + 'hcenter/hcenter_004.h5', 'r')['Data-Set-2'][:]
    # mpv.wplus[0][...] = h5py.File(base_filename + 'wplusx/wplusx_004.h5', 'r')['Data-Set-2'][:]
    # mpv.wplus[1][...] = h5py.File(base_filename + 'wplusy/wplusy_004.h5', 'r')['Data-Set-2'][:]
    if writer != None:
        writer.populate('after_ebnaimp','hcenter',mpv.wcenter)
        writer.populate('after_ebnaimp','wplusx',mpv.wplus[0])
        writer.populate('after_ebnaimp','wplusy',mpv.wplus[1])

    rhs[...], rhs_max = divergence_nodes(rhs,elem,node,Sol,ud)
    rhs /= dt
    
    if ud.is_compressible:
        rhs = rhs_from_p_old(rhs,node,mpv)

    x_wall = ud.bdry_type[0] == BdryType.WALL
    y_wall = ud.bdry_type[1] == BdryType.WALL

    if y_wall:
        rhs[:,2] *= 2.
        rhs[:,-3] *= 2.
    if x_wall:
        rhs[2,:] *= 2.
        rhs[-3,:] *= 2.

    # rhs = h5py.File(base_filename + 'rhs_nodes/rhs_nodes_004.h5','r')['Data-Set-2']
    if writer != None:
        writer.populate('after_ebnaimp','rhs_nodes',rhs)

    # lap2D = stencil_9pt_2nd_try(rhs,elem,node,mpv,ud)
    lap2D = stencil_9pt_3rd_try(elem,node,mpv,ud)

    sh = (ud.inx)*(ud.iny)

    lap2D = LinearOperator((sh,sh),matvec=lap2D)
    # lap_test = np.zeros_like(mpv.wcenter)
    # lap_test[node.igx:-node.igx,node.igy:-node.igy] = lap2D(p2.reshape(-1,)).reshape(257,22)

    # if writer != None:
    #     writer.populate('after_ebnaimp','lap_test',lap_test)

    p2,_ = bicgstab(lap2D,rhs[node.igx:-node.igx,node.igy:-node.igy].reshape(-1,),x0=p2.reshape(-1,),tol=1e-6,maxiter=500)

    p2_full = np.zeros(nc).squeeze()
    p2_full[node.igx:-node.igx,node.igy:-node.igy] = p2.reshape(ud.inx,ud.iny)

    # p2_full = h5py.File(base_filename + 'pnew/p2_full_004.h5', 'r')['Data-Set-2'][:]
    if writer != None:
        print("Yes?!")
        writer.populate('after_ebnaimp','p2_full',p2_full)

    mpv.dp2_nodes[...] = np.copy(p2_full)
    correction_nodes(Sol,elem,node,mpv,p2_full,dt,ud)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    mpv.p2_nodes[...] = p2_full
    mpv.rhs = rhs
    set_ghostnodes_p2(mpv.p2_nodes,node,ud)


def correction_nodes(Sol,elem,node,mpv,p,dt,ud):
    ndim = node.ndim
    # coriolis = ud.coriolis_strength[0]
    # time_offset = 3.0 - ud.acoustic_order
    vcorr = 1.

    igx,igy,igz = node.igs[0],node.igs[1],node.igs[2]
    dx = node.dx
    dy = node.dy
    oodx = 1.0/dx
    oody = 1.0/dy

    hplusx = mpv.wplus[0]
    hplusy = mpv.wplus[1]

    # dSdy = (mpv.HydroState_n.S0[0,igy+1:-igy] - mpv.HydroState_n.S0[0,igy:-igy-1]) * oody
    inner_idx = (slice(igx,-igx),slice(igy,-igy))
    # inner_idx = (slice(1,-1),slice(1,-1))
    # inner_idx = (slice(None,None), slice(None,None))
    inner_eidx = (slice(igx,-igx-1),slice(igy,-igy-1))
    # inner_eeidx = (slice(igx+1,-igx-1),slice(igy+1,-igy-1))
    
    p = p[inner_idx]

    indices = [idx for idx in product([slice(0,-1),slice(1,None)], repeat=ndim)]
    print(indices)
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
    # print(Dpx.shape)
    # print(hplusx[inner_eeidx].shape)
    # print(Sol.rhou[inner_eidx].shape)

    # kernel_Dpx = np.array([[-1.,1],[-1.,1.]])
    # kernel_Dpy = np.array([[-1.,-1],[1.,1.]])

    # Dpx = 0.5 * oodx * signal.convolve2d(p,kernel_Dpx,mode='valid')
    # Dpy = 0.5 * oody * signal.convolve2d(p,kernel_Dpy,mode='valid')
    
    # find_nearest(Dpx,62915.298150945368)
    thinv = Sol.rho[inner_idx] / Sol.rhoY[inner_idx]

    Sol.rhou[inner_idx] += -dt * thinv * hplusx[inner_eidx] * Dpx
    Sol.rhov[inner_idx] += -dt * thinv * hplusy[inner_eidx] * Dpy * vcorr
    # Sol.rhoX[inner_idx] += -time_offset * dt * dSdy * Sol.rhov[inner_idx] * vcorr


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

    # idx = 521
    # idx_s = 0

    # x_periodic = ud.bdry_type[0] == BdryType.PERIODIC
    # y_periodic = ud.bdry_type[1] == BdryType.PERIODIC

    # icxe = elem.icx
    # icye = elem.icy

    # hplusx = np.zeros_like(mpv.wplus[0].reshape(-1,))
    # hplusy = np.zeros_like(mpv.wplus[1].reshape(-1,))

    # shp = mpv.wplus[0].shape
    # print(shp)

    # for j in range(igy - y_periodic, icye - igy + y_periodic):
    #     m = j * icxe
    #     strat = 2.0 * (mpv.HydroState_n.Y0.ravel()[j+1] - mpv.HydroState_n.Y0.ravel()[j]) / (mpv.HydroState_n.Y0.ravel()[j+1] + mpv.HydroState_n.Y0.ravel()[j]) / dy

    #     for i in range(igx - x_periodic, icxe - igx + x_periodic):
    #         n = m + i

    #         Y = Sol.rhoY.T.ravel()[n] / Sol.rho.T.ravel()[n]
    #         coeff = Gammainv * Sol.rhoY.T.ravel()[n] * Y
    #         fsqsc = dt * dt * coriolis * coriolis
    #         fimp = 1.0 / (1.0 + fsqsc)
    #         Nsqsc = time_offset * dt*dt * (g/Msq) * strat
    #         gimp = 1.0 / (nonhydro + Nsqsc)

    #         hplusx[n] = coeff * gimp
    #         hplusy[n] = coeff * fimp

    # mpv.wplus[0][...] = hplusx.reshape(26,261).T
    # mpv.wplus[1][...] = hplusy.reshape(26,261).T

    # coeff = coeff.reshape(-1,1)
    # print(Sol.rhoY.ravel()[14+13])
    # print(Y.ravel()[idx], coeff.ravel()[idx], fsqsc, fimp, Nsqsc.ravel()[idx_s], gimp.ravel()[idx_s], Sol.rhoY.ravel()[idx])
    # find_nearest(Sol.rhoY,0.9893096665665001)

    kernel = (ndim - 1) * np.ones((2,2))
    # kernel = np.ones((2,2))

    for iz in range(icz):
        iz = slice(None,) if iz == 0 else iz

        mpv.wcenter[iz][innerdim] = ccenter * signal.convolve2d(Sol.rhoY[igx-1:-igx+1,igy-1:-igy+1]**cexp,kernel,mode='valid') / kernel.sum()
    # mpv.wcenter[innerdim] = ccenter * signal.convolve2d(Sol.rhoY[igx-1:-igx+1,igy-1:-igy+1]**cexp,kernel,mode='valid') / kernel.sum()
        
    scale_wall_node_values(mpv.wcenter, node, ud)


def scale_wall_node_values(rhs, node, ud, factor=.5):
    # Test: should be correct.
    ndim = node.ndim
    igs = node.igs

    wall_idx = np.empty((ndim), dtype=object)
    for dim in range(ndim):
        wall_idx[dim] = slice(igs[dim],-igs[dim])

    for dim in range(ndim):
        is_wall = ud.bdry_type[dim] == BdryType.WALL
        if is_wall:
            for direction in [-1,1]:
                wall_idx[dim] = (igs[dim]) * direction
                if direction == -1:
                    wall_idx[dim] -= 1
                wall_idx_tuple = tuple(wall_idx)
                rhs[wall_idx_tuple] *= factor
            
            


def divergence_nodes(rhs,elem,node,Sol,ud):
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

    rhs_max = np.max(rhs[inner_idx]) if np.max(rhs[inner_idx]) > 0 else 0
    return rhs, rhs_max


def rhs_from_p_old(rhs,node,mpv):
    igs = node.igs
    ndim = node.ndim

    assert ndim != 1, "Not implemented for 1D"
    
    inner_idx = np.empty((ndim), dtype=object)
    for dim in range(ndim):
        inner_idx[dim] = slice(igs[dim],-igs[dim])
    
    inner_idx = tuple(inner_idx)
    rhs_hh = mpv.wcenter[inner_idx] * mpv.p2_nodes[inner_idx]
    rhs[inner_idx] += rhs_hh
    return rhs


