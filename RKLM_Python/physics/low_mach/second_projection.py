from inputs.enum_bdry import BdryType
from inputs.boundary import set_explicit_boundary_data, set_ghostnodes_p2
from physics.low_mach.laplacian import stencil_9pt, stencil_5pt, stencil_3pt, precon_diag_prepare
from scipy import signal
import numpy as np
from itertools import product

from scipy.sparse.linalg import LinearOperator, spsolve, bicgstab, gmres
from scipy.sparse import eye, diags
import matplotlib.pyplot as plt

from management.debug import find_nearest
import h5py

# taken from https://stackoverflow.com/questions/33512081/getting-the-number-of-iterations-of-scipys-gmres-iterative-method
class solver_counter(object):
    def __init__(self, disp=True):
        self.niter = 0
    def __call__(self, rk=None):
        self.niter += 1
        self.rk = rk
        # self.rk = rk[0]


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

    scale_wall_node_values(mpv.rhs, node, ud, 2.0)

    ## 2D-case ###
    inner_idx = (slice(1,-1),slice(1,-1))
    p2n = p2n[inner_idx]

    dpdx_kernel = np.array([[-1.,1.],[-1.,1.]])
    dpdy_kernel = np.array([[-1.,-1.],[1.,1.]])

    dpdx = -wp * 0.5 * signal.convolve2d(p2n, dpdx_kernel, mode='valid') / dx
    dpdy = -wp * 0.5 * signal.convolve2d(p2n, dpdy_kernel, mode='valid') / dy

    dpdx = dpdx.T
    dpdy = dpdy.T

    Sol.flip()

    igx = elem.igx
    igy = elem.igy

    dSdy = mpv.HydroState_n.S0[igy-1:-igy+1]
    dSdy = signal.convolve(dSdy,[1.,-1.],mode='valid').reshape(-1,1)
    dSdy = np.repeat(dSdy,elem.icx-igx,axis=1) / dx

    S0c = mpv.HydroState.S0[igy-1:-igy+1].reshape(-1,1)
    S0c = np.repeat(S0c,elem.icx-igx,axis=1)

    v = Sol.rhou / Sol.rho
    v = v[inner_idx]
    time_offset_expl = ud.acoustic_order - 1.0

    rhoYovG = Ginv * Sol.rhoY[inner_idx]
    dchi = Sol.rhoX[inner_idx] / Sol.rho[inner_idx]
    dbuoy = -1. * Sol.rhoY[inner_idx] * dchi
    drhou = Sol.rhov[inner_idx] - u0 * Sol.rho[inner_idx]

    rhoY = Sol.rhoY**(th.gamm - 2.0)
    dpidP_kernel = np.array([[1.,1.],[1.,1.]])
    dpidP = (th.gm1 / ud.Msq) * signal.convolve2d(rhoY, dpidP_kernel, mode='valid') / dpidP_kernel.sum()

    Sol.rhov[inner_idx] = Sol.rhov[inner_idx] + dt * ( -1. * rhoYovG * dpdy + coriolis * Sol.rhow[inner_idx])

    Sol.rhou[inner_idx] = Sol.rhou[inner_idx] + dt * ( -1. * rhoYovG * dpdx + (g/Msq) * dbuoy) * nonhydro * (1 - ud.is_ArakawaKonor)

    Sol.rhow[inner_idx] = Sol.rhow[inner_idx] - dt * coriolis * drhou

    Sol.rhoX[inner_idx] = (Sol.rho[inner_idx] * (Sol.rho[inner_idx] / Sol.rhoY[inner_idx] - S0c)) + time_offset_expl * dt * (-v * dSdy) * Sol.rho[inner_idx]

    dp2n[inner_idx] -= dt * dpidP * div[inner_idx]

    if (ud.is_compressible != 0):
        weight = ud.acoustic_order - 1.0
        mpv.p2_nodes += weight * dp2n.T

    Sol.flip()

    set_ghostnodes_p2(mpv.p2_nodes,node, ud)
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

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

    # first_nodes_row_right_idx = (slice(0,1), slice(1,None))
    # first_nodes_row_left_idx = (slice(0,1),slice(0,-1))

    first_nodes_row_right_idx = (slice(1,None))
    first_nodes_row_left_idx = (slice(0,-1))

    strat = 2.0 * (mpv.HydroState_n.Y0[first_nodes_row_right_idx] - mpv.HydroState_n.Y0[first_nodes_row_left_idx]) / (mpv.HydroState_n.Y0[first_nodes_row_right_idx] + mpv.HydroState_n.Y0[first_nodes_row_left_idx])

    strat = strat.reshape(1,-1)
    strat = np.repeat(strat, elem.icx, axis=0) / dy

    Nsqsc = time_offset * dt**2 * (g / Msq) * strat
    
    dbuoy = Sol.rhoY * Sol.rhoX / Sol.rho
    rhov = (nonhydro * Sol.rhov - dt * (g/Msq) * dbuoy) / (nonhydro + Nsqsc)

    drhou = Sol.rhou - u0 * Sol.rho
    Sol.rhou[...] = u0 * Sol.rho + ooopfsqsc * (drhou + dt * coriolis * Sol.rhow)
    Sol.rhov[...] = rhov
    Sol.rhow[...] = ooopfsqsc * (Sol.rhow - dt * coriolis * drhou)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)


total_iter = 0
total_calls = 0
def euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, dt, alpha_diff, writer = None, label=None):
    nc = node.sc
    rhs = np.zeros_like(mpv.p2_nodes)

    p2 = np.copy(mpv.p2_nodes[node.igx:-node.igx,node.igy:-node.igy])
    
    if writer != None:
        writer.populate(str(label)+'_after_ebnaimp','p2_initial',mpv.p2_nodes)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt)
    
    if writer != None:
        writer.populate(str(label)+'_after_ebnaimp','hcenter',mpv.wcenter)
        writer.populate(str(label)+'_after_ebnaimp','wplusx',mpv.wplus[0])
        writer.populate(str(label)+'_after_ebnaimp','wplusy',mpv.wplus[1])

    rhs[...], _ = divergence_nodes(rhs,elem,node,Sol,ud)
    rhs /= dt

    if ud.is_compressible == 1:
        rhs = rhs_from_p_old(rhs,node,mpv)
    elif ud.is_compressible == 0:
        if ud.is_ArakawaKonor:
            rhs -= mpv.wcenter * mpv.dp2_nodes
            mpv.wcenter[...] = 0.0
        else:
            mpv.wcenter[...] = 0.0
    # else:
        # mpv.wcenter *= ud.compressibility

    diag_inv = precon_diag_prepare(mpv, elem, node, ud)
    rhs *= diag_inv

    # if y_wall:
    #     rhs[:,2] *= 2.
    #     rhs[:,-3] *= 2.
    # if x_wall:
    #     rhs[2,:] *= 2.
    #     rhs[-3,:] *= 2.

    if writer != None:
        writer.populate(str(label)+'_after_ebnaimp','rhs_nodes',rhs)

    lap2D = stencil_9pt(elem,node,mpv,ud,diag_inv)

    sh = (ud.inx)*(ud.iny)

    lap2D = LinearOperator((sh,sh),lap2D)
    
    counter = solver_counter()
    
    p2,info = bicgstab(lap2D,rhs[node.igx:-node.igx,node.igy:-node.igy].ravel(),x0=p2.ravel(),tol=1e-10,maxiter=6000,callback=counter)

    # print("Convergence info = %i, no. of iterations = %i" %(info,counter.niter))

    global total_calls, total_iter
    total_iter += counter.niter
    total_calls += 1
    # print("Total calls to BiCGStab routine = %i, total iterations = %i" %(total_calls, total_iter))

    p2_full = np.zeros(nc).squeeze()
    p2_full[node.igx:-node.igx,node.igy:-node.igy] = p2.reshape(ud.inx,ud.iny)

    if writer != None:
        writer.populate(str(label)+'_after_ebnaimp','p2_full',p2_full)

    correction_nodes(Sol,elem,node,mpv,p2_full,dt,ud)
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    mpv.dp2_nodes[...] = p2_full - mpv.p2_nodes
    
    # if ud.is_ArakawaKonor:
    #     if ud.is_compressible == 0:
    #         dp2_rhs = exner_perturbation_constraint(Sol,elem,th,p2.reshape(ud.inx,ud.iny))

    #         # lap2D_exner = stencil_3pt(elem,node,ud)
    #         lap2D_exner = stencil_5pt(elem,node,ud)
    #         sh = (ud.inx)*(ud.iny)

    #         lap2D_exner = LinearOperator((sh,sh),lap2D_exner)

    #         dp2,_ = bicgstab(lap2D_exner,dp2_rhs.ravel(),tol=1e-16,maxiter=6000)

    #         mpv.dp2_nodes[node.igx:-node.igx,node.igy:-node.igy] = dp2.reshape(ud.inx,ud.iny)
            
    mpv.p2_nodes[...] = p2_full
    mpv.rhs = rhs
    set_ghostnodes_p2(mpv.dp2_nodes,node,ud)
    set_ghostnodes_p2(mpv.p2_nodes,node,ud)


def exner_perturbation_constraint(Sol,elem,th,p2):
    # 2D function!
    rhs = np.zeros_like(Sol.rhou)

    gathering_kernel = np.array([[1.,1.],[1.,1.]])
    gathering_kernel /= gathering_kernel.sum()
    rhs = signal.convolve2d(rhs, gathering_kernel, mode='valid')

    return rhs[1:-1,1:-1]


def correction_nodes(Sol,elem,node,mpv,p,dt,ud):
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

    dSdy = (mpv.HydroState_n.S0[igy+1:-igy] - mpv.HydroState_n.S0[igy:-igy-1]) * oody
    inner_idx = (slice(igx,-igx),slice(igy,-igy))
    inner_eidx = (slice(igx,-igx-1),slice(igy,-igy-1))
    
    p = p[inner_idx]

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

    # ccenter = - (ud.compressibility * ud.Msq) * th.gm1inv / (dt**2) / time_offset
    ccenter = - ud.Msq * th.gm1inv / (dt**2) / time_offset
    cexp = 2.0 - th.gamm

    igs = elem.igs
    igx = igs[0]
    igy = igs[1]

    nindim = np.empty((ndim),dtype='object')
    innerdim = np.copy(nindim)
    eindim = np.empty((ndim),dtype='object')

    for dim in range(ndim):
        is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
        nindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
        innerdim[dim] = slice(igs[dim],-igs[dim])
        eindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)

        if dim == 1:
            y_idx = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)
            right_idx = None if -igs[dim]+is_periodic == 0 else -igs[dim]+is_periodic
            y_idx1 = slice(igs[dim]-is_periodic+1, right_idx)
 
    strat = 2.0 * (mpv.HydroState_n.Y0[y_idx1] - mpv.HydroState_n.Y0[y_idx]) / (mpv.HydroState_n.Y0[y_idx1] + mpv.HydroState_n.Y0[y_idx]) / dy

    nindim = tuple(nindim)
    eindim = tuple(eindim)
    innerdim = tuple(innerdim)

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

    kernel = (ndim - 1) * np.ones((2,2))

    for iz in range(icz):
        iz = slice(None,) if iz == 0 else iz

        mpv.wcenter[iz][innerdim] = ccenter * signal.convolve2d(Sol.rhoY[igx-1:-igx+1,igy-1:-igy+1]**cexp,kernel,mode='valid') / kernel.sum()
        
    scale_wall_node_values(mpv.wcenter, node, ud)


def scale_wall_node_values(rhs, node, ud, factor=.5):
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
