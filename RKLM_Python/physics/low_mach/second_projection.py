from inputs.enum_bdry import BdryType
from inputs.boundary import set_explicit_boundary_data, set_ghostnodes_p2
from physics.low_mach.laplacian import stencil_9pt, stencil_32pt, stencil_hs, precon_diag_prepare
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

def euler_forward_non_advective(Sol, mpv, elem, node, dt, ud, th):
    nonhydro = ud.nonhydrostasy
    wp = 1. # with pressure

    g = ud.gravity_strength[1]
    Msq = ud.Msq
    Ginv = th.Gammainv
    coriolis = ud.coriolis_strength[0]
    u0 = ud.wind_speed

    p2n = np.copy(mpv.p2_nodes)
    dp2n = np.zeros_like(p2n)
    ndim = elem.ndim
    dx, dy, dz = node.dx, node.dy, node.dz

    mpv.rhs = np.array(mpv.rhs)
    mpv.rhs *= 0.0
    mpv.rhs[...], _ = divergence_nodes(mpv.rhs,elem,node,Sol,ud)
    # div = mpv.rhs

    scale_wall_node_values(mpv.rhs, node, ud, 2.0)
    div = mpv.rhs

    i1 = np.empty(ndim, dtype='object')
    for dim in range(ndim):
        i1[dim] = slice(1,-1)
    i1 = tuple(i1)
    p2n = p2n[i1]

    # kernels = np.empty(ndim, dtype='ndarray')
    kernels = []
    for dim in range(ndim):
        kernel = np.ones([2]*ndim)
        slc = [slice(None,)]*ndim
        slc[dim] = slice(0,1)
        kernel[tuple(slc)] *= -1.0
        kernels.append(kernel)

    # dpdy_kernel = np.array([[-1.,1.],[-1.,1.]])
    # dpdx_kernel = np.array([[-1.,-1.],[1.,1.]])

    # dpdx = -wp * 0.5 * signal.convolve2d(p2n, dpdx_kernel, mode='valid') / dx
    # dpdy = -wp * 0.5 * signal.convolve2d(p2n, dpdy_kernel, mode='valid') / dy

    dpdx = -wp * 0.5**(ndim-1) * signal.fftconvolve(p2n, kernels[0], mode='valid') / dx
    dpdy = -wp * 0.5**(ndim-1) * signal.fftconvolve(p2n, kernels[1], mode='valid') / dy
    if ndim == 3: dpdz = -wp * 0.5**(ndim-1) * signal.fftconvolve(p2n, kernels[2], mode='valid') / dz
    else: dpdz = 0

    igx, igy, igz, igs = elem.igx, elem.igy, elem.igz, elem.igs

    dSdy = mpv.HydroState_n.S0[igy-1:-igy+1]
    dSdy = signal.convolve(dSdy,[1.,-1.],mode='valid') / dy
    # dSdy = np.repeat(dSdy,elem.icx-igx,axis=1).T / dy

    S0c = mpv.HydroState.S0[igy-1:-igy+1]

    for dim in range(0,ndim,2):
        dSdy = np.expand_dims(dSdy, dim)
        dSdy = np.repeat(dSdy, elem.sc[dim]-igs[dim], axis=dim)
        S0c = np.expand_dims(S0c, dim)
        S0c = np.repeat(S0c, elem.sc[dim]-igs[dim], axis=dim)

    v = Sol.rhov / Sol.rho
    v = v[i1]
    time_offset_expl = ud.acoustic_order - 1.0

    rhoYovG = Ginv * Sol.rhoY[i1]
    dchi = Sol.rhoX[i1] / Sol.rho[i1]
    dbuoy = -1. * Sol.rhoY[i1] * dchi
    drhou = Sol.rhou[i1] - u0 * Sol.rho[i1]

    rhoY = Sol.rhoY**(th.gamm - 2.0)
    dpidP_kernel = np.ones([2]*ndim)
    dpidP = (th.gm1 / ud.Msq) * signal.fftconvolve(rhoY, dpidP_kernel, mode='valid') / dpidP_kernel.sum()

    Sol.rhou[i1] = Sol.rhou[i1] + dt * ( -1. * rhoYovG * dpdx + coriolis * Sol.rhow[i1])

    Sol.rhov[i1] = Sol.rhov[i1] + dt * ( -1. * rhoYovG * dpdy + (g/Msq) * dbuoy) * nonhydro * (1 - ud.is_ArakawaKonor)

    Sol.rhow[i1] = Sol.rhow[i1] + dt * ( -(ndim == 3) * rhoYovG * dpdz - coriolis * drhou)

    Sol.rhoX[i1] = (Sol.rho[i1] * (Sol.rho[i1] / Sol.rhoY[i1] - S0c)) + time_offset_expl * dt * (-v * dSdy) * Sol.rho[i1]

    dp2n[i1] -= dt * dpidP * div[i1]

    # if (ud.is_compressible == 1):
    weight = ud.compressibility * ud.acoustic_order - 1.0
    mpv.p2_nodes += weight * dp2n

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

    first_nodes_row_right_idx = (slice(1,None))
    first_nodes_row_left_idx = (slice(0,-1))

    strat = 2.0 * (mpv.HydroState_n.Y0[first_nodes_row_right_idx] - mpv.HydroState_n.Y0[first_nodes_row_left_idx]) / (mpv.HydroState_n.Y0[first_nodes_row_right_idx] + mpv.HydroState_n.Y0[first_nodes_row_left_idx])

    s0 = mpv.HydroState.S0

    for dim in range(0,elem.ndim,2):
        strat = np.expand_dims(strat, dim)
        strat = np.repeat(strat, elem.sc[dim], axis=dim)
        s0 = np.expand_dims(s0, dim)
        s0 = np.repeat(s0, elem.sc[dim], axis=dim)
    strat /= dy

    Nsqsc = time_offset * dt**2 * (g / Msq) * strat

    dbuoy = -Sol.rhoY * (Sol.rhoX / Sol.rho)
    rhov = (nonhydro * Sol.rhov + dt * (g/Msq) * dbuoy) / (nonhydro + Nsqsc)

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

    if elem.ndim == 2:
        p2 = np.copy(mpv.p2_nodes[node.igx:-node.igx,node.igy:-node.igy])
    elif elem.ndim == 3:
        # p2 = np.copy(mpv.p2_nodes[node.igx:-node.igx,node.igy:-node.igy,node.igz:-node.igz])
        p2 = np.copy(mpv.p2_nodes[1:-1,1:-1,1:-1])

        # p2 = np.copy(mpv.p2_nodes)
    
    if writer != None:
        writer.populate(str(label),'p2_initial',mpv.p2_nodes)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt)

    i0 = node.ndim * [(slice(0,-1))]
    i0 = tuple(i0)

    if writer != None:
        writer.populate(str(label),'hcenter',mpv.wcenter)
        writer.populate(str(label),'wplusx',mpv.wplus[0])
        writer.populate(str(label),'wplusy',mpv.wplus[1])
        if elem.ndim == 3: writer.populate(str(label),'wplusz',mpv.wplus[2])
        else: writer.populate(str(label),'wplusz',np.zeros_like(mpv.wplus[0]))

    rhs[...], _ = divergence_nodes(rhs,elem,node,Sol,ud)
    rhs /= dt

    if writer != None:
        writer.populate(str(label),'rhs',rhs)

    if ud.is_compressible == 1:
        rhs = rhs_from_p_old(rhs,node,mpv)
    # if 
    elif ud.is_compressible == 0:
        if ud.is_ArakawaKonor:
            rhs -= mpv.wcenter * mpv.dp2_nodes
            mpv.wcenter[...] = 0.0
        else:
            rhs_new = rhs_from_p_old(rhs,node,mpv)
            rhs = ud.compressibility * rhs_new + (1.0 - ud.compressibility) * rhs
            mpv.wcenter[...] *= ud.compressibility
    else:
        mpv.wcenter *= ud.compressibility

    if writer != None:
        writer.populate(str(label),'rhs_nodes',rhs)

    mpv.rhs = rhs

    diag_inv = precon_diag_prepare(mpv, elem, node, ud)
    rhs *= diag_inv

    if elem.ndim == 2:
        lap = stencil_9pt(elem,node,mpv,ud,diag_inv)

        sh = (ud.inx)*(ud.iny)

        # lap = LinearOperator((sh,sh),lap)
    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] > 2:
        lap = stencil_32pt(elem,node,mpv,ud,diag_inv,dt)
        # p2 = mpv.p2_nodes[1:-1,1:-1,1:-1]

        sh = p2.reshape(-1).shape[0]
    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] <= 2:
        p2 = np.copy(mpv.p2_nodes[1:-1,elem.igs[1],1:-1])
        lap = stencil_hs(elem,node,mpv,ud,diag_inv,dt)
        sh = p2.reshape(-1).shape[0]
    lap = LinearOperator((sh,sh),lap)
    
    counter = solver_counter()

    if elem.ndim == 2:
        rhs_inner = rhs[node.igx:-node.igx,node.igy:-node.igy].ravel()
    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] > 2:
        rhs_inner = rhs[1:-1,1:-1,1:-1].ravel()
        # rhs_inner = rhs.ravel()
        # p2 = mpv.p2_nodes[1:-1,1:-1,1:-1]
        # rhs_inner = rhs.ravel()
    else:
        rhs_inner = rhs[1:-1,elem.igs[1],1:-1].ravel()
    p2,info = bicgstab(lap,rhs_inner,tol=ud.tol,maxiter=ud.max_iterations,callback=counter)
    # p2,info = bicgstab(lap,rhs.ravel(),x0=p2.ravel(),tol=1e-16,maxiter=6000,callback=counter)
    # print("Convergence info = %i, no. of iterations = %i" %(info,counter.niter))

    global total_calls, total_iter
    total_iter += counter.niter
    total_calls += 1
    # print(counter.niter)
    # print("Total calls to BiCGStab routine = %i, total iterations = %i" %(total_calls, total_iter))

    p2_full = np.zeros(nc).squeeze()
    if elem.ndim == 2:
        p2_full[node.igx:-node.igx,node.igy:-node.igy] = p2.reshape(ud.inx,ud.iny)
    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] > 2:
        p2_full[1:-1,1:-1,1:-1] = p2.reshape(ud.inx+2,ud.iny+2,ud.inz+2)
    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] <= 2:
        p2 = p2.reshape(ud.inx+2, ud.inz+2)
        p2 = np.expand_dims(p2,1)
        p2 = np.repeat(p2, node.icy, axis=1)
        p2_full[1:-1,:,1:-1] = p2
    if writer != None:
        writer.populate(str(label),'p2_full',p2_full)

    correction_nodes(Sol,elem,node,mpv,p2_full,dt,ud)
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    mpv.p2_nodes[...] = p2_full
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
    coriolis = ud.coriolis_strength[0]
    time_offset = 3.0 - ud.acoustic_order

    igs, igy = node.igs, node.igy
    oodxyz = 1.0 / node.dxyz
    oodx , oody ,oodz = oodxyz[0], oodxyz[1], oodxyz[2]

    hplusx = mpv.wplus[0]
    hplusy = mpv.wplus[1]
    if ndim == 3: hplusz = mpv.wplus[2]

    dSdy = (mpv.HydroState_n.S0[igy+1:-igy] - mpv.HydroState_n.S0[igy:-igy-1]) * oody

    for dim in range(0,ndim,2):
        dSdy = np.expand_dims(dSdy, dim)
        dSdy = np.repeat(dSdy, elem.sc[dim]-2*igs[dim], axis=dim)

    n2e, i2 = np.empty(ndim, dtype='object'), np.empty(ndim, dtype='object')
    for dim in range(ndim):
        n2e[dim] = slice(0,-1)
        i2[dim] = slice(igs[dim],-igs[dim])
    n2e, i2 = tuple(n2e), tuple(i2)
    
    p = p[i2]

    indices = [idx for idx in product([slice(0,-1),slice(1,None)], repeat=ndim)]
    if ndim == 2:
        signs_x = (-1., -1., +1., +1.)
        signs_y = (-1., +1., -1., +1.)
        signs_z = ( 0.,  0.,  0.,  0.)
    elif ndim == 3:
        signs_x = (-1., -1., -1., -1., +1., +1., +1., +1.)
        signs_y = (-1., -1., +1., +1., -1., -1., +1., +1.)
        signs_z = (-1., +1., -1., +1., -1., +1., -1., +1.)

    Dpx, Dpy, Dpz = 0., 0., 0.
    cnt = 0
    for index in indices:
        Dpx += signs_x[cnt] * p[index]
        Dpy += signs_y[cnt] * p[index]
        Dpz += signs_z[cnt] * p[index]
        cnt += 1

    Dpx *= 0.5**(ndim-1) * oodx
    Dpy *= 0.5**(ndim-1) * oody
    Dpz *= 0.5**(ndim-1) * oodz

    thinv = Sol.rho[i2] / Sol.rhoY[i2]

    Sol.rhou[i2] += -dt * thinv * hplusx[n2e][i2] * (Dpx + dt * coriolis * Dpz)
    Sol.rhov[i2] += -dt * thinv * hplusy[n2e][i2] * Dpy
    if ndim == 3: Sol.rhow[i2] += -dt * thinv * hplusz[n2e][i2] * (Dpz - dt * coriolis * Dpx)
    Sol.rhoX[i2] += -time_offset * dt * dSdy * Sol.rhov[i2]


def operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt):
    g = ud.gravity_strength[1]
    Msq = ud.Msq
    Gammainv = th.Gammainv

    ndim = node.ndim
    nonhydro = ud.nonhydrostasy
    dy = elem.dy

    time_offset = 3.0 - ud.acoustic_order

    coriolis = ud.coriolis_strength[0]

    ccenter = - ud.Msq * th.gm1inv / (dt**2) / time_offset
    cexp = 2.0 - th.gamm

    igs = elem.igs
    igx = igs[0]
    igy = igs[1]

    nindim = np.empty((ndim),dtype='object')
    innerdim = np.copy(nindim)
    innerdim1 = np.copy(nindim)
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

        innerdim1[dim] = slice(igs[dim]-1, (-igs[dim]+1))
 
    strat = 2.0 * (mpv.HydroState_n.Y0[y_idx1] - mpv.HydroState_n.Y0[y_idx]) / (mpv.HydroState_n.Y0[y_idx1] + mpv.HydroState_n.Y0[y_idx]) / dy
    # strat = 2.0 * (np.diff(mpv.HydroState_n.Y0)[igs[1]:-igs[1]]) / np.diff(mpv.HydroState_n.Y0)[igs[1]:-igs[1]] / dy

    nindim = tuple(nindim)
    eindim = tuple(eindim)
    innerdim = tuple(innerdim)
    innerdim1 = tuple(innerdim1)

    for dim in range(0,elem.ndim,2):
        is_periodic = ud.bdry_type[dim] != BdryType.PERIODIC
        strat = np.expand_dims(strat, dim)
        strat = np.repeat(strat, elem.sc[dim]-int(2*is_periodic+igs[dim]), axis=dim)

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

    kernel = np.ones([2] * ndim)

    mpv.wcenter[innerdim] = ccenter * signal.fftconvolve(Sol.rhoY[innerdim1]**cexp,kernel,mode='valid') / kernel.sum()

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
    dxyz = node.dxyz
    inner_idx = np.empty((ndim), dtype=object)

    for dim in range(ndim):
        is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
        inner_idx[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
        
    indices = [idx for idx in product([slice(0,-1),slice(1,None)], repeat = ndim)]
    signs = [sgn for sgn in product([1,-1], repeat = ndim)]
    inner_idx = tuple(inner_idx)
    Sols = np.stack((Sol.rhou[inner_idx], Sol.rhov[inner_idx], Sol.rhow[inner_idx]), axis=-1)

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
    rhs_n = np.zeros_like(rhs)
    rhs_hh = mpv.wcenter[inner_idx] * mpv.p2_nodes[inner_idx]
    rhs_n[inner_idx] = rhs[inner_idx] + rhs_hh
    return rhs_n
