from inputs.enum_bdry import BdryType
from inputs.boundary import set_explicit_boundary_data, set_ghostnodes_p2
from physics.low_mach.laplacian import stencil_9pt, stencil_27pt, stencil_hs, stencil_vs, stencil_9pt_numba_test, precon_diag_prepare
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

def euler_forward_non_advective(Sol, mpv, elem, node, dt, ud, th, writer = None, label = None):
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[1]
    Msq = ud.Msq
    Ginv = th.Gammainv
    corr_h1 = ud.coriolis_strength[0]
    corr_v  = ud.coriolis_strength[1]
    corr_h2 = ud.coriolis_strength[2]
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed

    p2n = np.copy(mpv.p2_nodes)
    dp2n = np.zeros_like(p2n)
    ndim = elem.ndim

    S0c = mpv.HydroState.get_S0c(elem)
    dSdy = mpv.HydroState_n.get_dSdy(elem,node)

    mpv.rhs[...] = 0.0
    mpv.rhs[...], _ = divergence_nodes(mpv.rhs,elem,node,Sol,ud)
    scale_wall_node_values(mpv.rhs, node, ud, 2.0)
    div = mpv.rhs

    if writer != None: writer.populate(str(label), 'rhs', div)

    v = Sol.rhov / Sol.rho

    rhoYovG = Ginv * Sol.rhoY
    dchi = Sol.rhoX / Sol.rho
    dbuoy = Sol.rhoY * dchi
    drhou = Sol.rhou - u0 * Sol.rho
    drhov = Sol.rhov - v0 * Sol.rho
    drhow = Sol.rhow - w0 * Sol.rho

    rhoY = Sol.rhoY**(th.gamm - 2.0)
    dpidP_kernel = np.ones([2]*ndim)
    dpidP = (th.gm1 / ud.Msq) * signal.fftconvolve(rhoY, dpidP_kernel, mode='valid') / dpidP_kernel.sum()

    dpdx, dpdy, dpdz = grad_nodes(p2n, elem, node)

    Sol.rhou = Sol.rhou - dt * ( rhoYovG * dpdx - corr_h2 * drhov + corr_v * drhow)

    Sol.rhov = Sol.rhov - dt * ( rhoYovG * dpdy + (g/Msq) * dbuoy * nonhydro - corr_h1 * drhow + corr_h2 * drhou) * (1 - ud.is_ArakawaKonor)

    Sol.rhow = Sol.rhow - dt * ( rhoYovG * dpdz - corr_v * drhou + corr_h1 * drhov) * (ndim == 3)

    Sol.rhoX = (Sol.rho * (Sol.rho / Sol.rhoY - S0c)) - dt * (v * dSdy) * Sol.rho

    dp2n[node.i1] -= dt * dpidP * div[node.i1]

    weight = ud.compressibility
    mpv.p2_nodes += weight * dp2n

    set_ghostnodes_p2(mpv.p2_nodes,node, ud)
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)


def euler_backward_non_advective_expl_part(Sol, mpv, elem, dt, ud, th):
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[1]
    Msq = ud.Msq

    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed

    wh1, wv, wh2 = dt * ud.coriolis_strength

    strat = mpv.HydroState_n.get_dSdy(elem, elem)

    Y = Sol.rhoY / Sol.rho
    nu = -dt**2 * (g / Msq) * strat * Y
    setattr(mpv,'nu',nu)

    drhou = Sol.rhou - u0 * Sol.rho
    drhov = Sol.rhov - v0 * Sol.rho
    drhow = Sol.rhow - w0 * Sol.rho

    # get coefficients of the explicit terms
    # common denominator
    denom = 1.0 / (wh1**2 + wh2**2 + (nu + nonhydro) * (wv**2 + 1.0))

    # U update
    coeff_uu = (wh1**2 + nu + nonhydro)
    coeff_uv = nonhydro * (wh1 * wv + wh2)
    coeff_uw = (wh1 * wh2 - (nu + nonhydro) * wv)

    # V update
    coeff_vu = (wh1 * wv - wh2)
    coeff_vv = nonhydro * (1 + wv**2)
    coeff_vw = (wh2 * wv + wh1)

    # W update
    coeff_wu = (wh1 * wh2 + (nu + nonhydro) * wv)
    coeff_wv = nonhydro * (wh2 * wv - wh1)
    coeff_ww = (nu + nonhydro + wh2**2)

    # Do the updates
    rhou = u0 * Sol.rho + denom * (coeff_uu * drhou + coeff_uv * drhov + coeff_uw * drhow)

    rhov = v0 * Sol.rho + denom * (coeff_vu * drhou + coeff_vv * drhov + coeff_vw * drhow)

    rhow = w0 * Sol.rho + denom * (coeff_wu * drhou + coeff_wv * drhov + coeff_ww * drhow)

    # Set the quantities
    Sol.rhou[...] = rhou
    Sol.rhov[...] = rhov
    Sol.rhow[...] = rhow

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)


total_iter = 0
total_calls = 0
def euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, dt, alpha_diff, Sol0 = None, writer = None, label=None):
    nc = node.sc
    rhs = np.zeros_like(mpv.p2_nodes)

    if elem.ndim == 2:
        p2 = np.copy(mpv.p2_nodes[node.igx:-node.igx,node.igy:-node.igy])
    elif elem.ndim == 3:
        p2 = np.copy(mpv.p2_nodes[1:-1,1:-1,1:-1])
    
    if writer != None:
        writer.populate(str(label),'p2_initial',mpv.p2_nodes)

    if Sol0 is not None:
        set_explicit_boundary_data(Sol0, elem, ud, th, mpv)
        operator_coefficients_nodes(elem, node, Sol0, mpv, ud, th, dt)
    else:
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

    set_ghostnodes_p2(mpv.p2_nodes,node,ud)
    correction_nodes(Sol,elem,node,mpv,mpv.p2_nodes,dt,ud,th,0)
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    rhs[...], _ = divergence_nodes(rhs,elem,node,Sol,ud)

    if writer != None:
        writer.populate(str(label),'rhs',rhs)

    rhs /= dt

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
    
    # diag_inv = np.ones_like(mpv.rhs)

    VS = True

    # prepare initial left-hand side and the laplacian stencil
    if elem.ndim == 2:
        lap = stencil_9pt(elem,node,mpv,Sol,ud,diag_inv,dt)

        # p2 = np.copy(mpv.p2_nodes[1:-1,2:-2])
        # shp = p2.shape
        # lap = stencil_9pt_numba_test(elem,node,mpv,Sol,ud,diag_inv,dt,th,shp)
        sh = (ud.inx)*(ud.iny)
        # p2 = np.copy(mpv.p2_nodes[1:-1,1:-1])
        # sh = p2.reshape(-1).shape[0]

    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] <= 2: 
        # horizontal slice hack
        p2 = np.copy(mpv.p2_nodes[1:-1,elem.igs[1],1:-1])
        lap = stencil_hs(elem,node,mpv,ud,diag_inv,dt)
        sh = p2.reshape(-1).shape[0]

    elif elem.ndim == 3 and elem.iicy > 1 and elem.iicz == 1:
        # vertical slice hack

        if not VS:
            p2 = np.copy(mpv.p2_nodes[...,elem.igz][node.igx:-node.igx,node.igy:-node.igy])
            lap = stencil_vs(elem,node,mpv,ud,diag_inv,dt)
            sh = (node.iicx)*(node.iicy)
        if VS:
            p2 = np.copy(mpv.p2_nodes[1:-1,1:-1,elem.igs[2]])
            lap = stencil_vs(elem,node,mpv,ud,diag_inv,dt)
            sh = p2.reshape(-1).shape[0]

    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] > 2:
        lap = stencil_27pt(elem,node,mpv,ud,diag_inv,dt)
        sh = p2.reshape(-1).shape[0]


    lap = LinearOperator((sh,sh),lap)
    # lap = LinearOperator(sh,lap)
    
    counter = solver_counter()

    # prepare right-hand side
    if elem.ndim == 2:
        rhs_inner = rhs[node.igx:-node.igx,node.igy:-node.igy].ravel()
        # rhs_inner = rhs[1:-1,1:-1].ravel()
    elif elem.ndim == 3 and elem.iicy > 1 and elem.iicz == 1:

        if not VS:
            rhs_inner = rhs[...,elem.igs[2]][node.igx:-node.igx,node.igy:-node.igy].ravel()
        if VS:
            rhs_inner = rhs[1:-1,1:-1,elem.igs[2]].ravel()

    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] > 2:
        rhs_inner = rhs[1:-1,1:-1,1:-1].ravel()
    else:
        rhs_inner = rhs[1:-1,elem.igs[1],1:-1].ravel()

    p2, _ = bicgstab(lap,rhs_inner,tol=ud.tol,maxiter=ud.max_iterations,callback=counter)
    # p2, _ = gmres(lap,rhs_inner,tol=ud.tol,maxiter=ud.max_iterations)
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
        # p2 = p2.reshape(ud.inx+2, ud.iny+2)
        # p2_full[1:-1,1:-1] = p2
    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] <= 2:
        # horizontal slice hack
        p2 = p2.reshape(ud.inx+2, ud.inz+2)
        p2 = np.expand_dims(p2,1)
        p2 = np.repeat(p2, node.icy, axis=1)
        p2_full[1:-1,:,1:-1] = p2
    elif elem.ndim == 3 and elem.iicy > 1 and elem.iicz == 1:
        if not VS:
            p2 = p2.reshape(node.iicx,node.iicy)
            p2 = np.repeat(p2[...,np.newaxis], node.icz, axis=2)
            p2_full[node.igx:-node.igx,node.igy:-node.igy] = p2
        if VS:
            p2 = p2.reshape(ud.inx+2, ud.iny+2)
            p2 = np.expand_dims(p2,2)
            p2 = np.repeat(p2, node.icz, axis=2)
            p2_full[1:-1,1:-1,:] = p2


    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] > 2:
        p2_full[1:-1,1:-1,1:-1] = p2.reshape(ud.inx+2,ud.iny+2,ud.inz+2)
    if writer != None:
        writer.populate(str(label),'p2_full',p2_full)

    set_ghostnodes_p2(p2_full,node,ud)
    correction_nodes(Sol,elem,node,mpv,p2_full,dt,ud,th,1)

    mpv.p2_nodes[...] += p2_full
    set_ghostnodes_p2(mpv.p2_nodes,node,ud)
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

 
    

def exner_perturbation_constraint(Sol,elem,th,p2):
    # 2D function!
    rhs = np.zeros_like(Sol.rhou)

    gathering_kernel = np.array([[1.,1.],[1.,1.]])
    gathering_kernel /= gathering_kernel.sum()
    rhs = signal.convolve2d(rhs, gathering_kernel, mode='valid')

    return rhs[1:-1,1:-1]


def correction_nodes(Sol,elem,node,mpv,p,dt,ud,th,updt_chi):
    ndim = node.ndim
    wh1, wv, wh2 = dt * ud.coriolis_strength
    nu = mpv.nu_c
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[1]
    Msq = ud.Msq

    Gammainv = th.Gammainv

    igs, igy = node.igs, node.igy
    oodxyz = 1.0 / node.dxyz
    oodx , oody ,oodz = oodxyz[0], oodxyz[1], oodxyz[2]

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
    thinv = 1.0
    coeff = Gammainv * Sol.rhoY[i2] #* Y
    Y = Sol.rhoY[i2] / Sol.rho[i2]

    nu = nu[n2e][i2]

    # get coefficients of the correction terms
    denom = 1.0 / (wh1**2 + wh2**2 + (nu + nonhydro) * (wv**2 + 1.0))
    
    # U update
    coeff_uu = (wh1**2 + nu + nonhydro)
    coeff_uv = nonhydro * (wh1 * wv + wh2)
    coeff_uw = (wh1 * wh2 - (nu + nonhydro) * wv)

    # V update
    coeff_vu = (wh1 * wv - wh2)
    coeff_vv = nonhydro * (1 + wv**2)
    coeff_vw = (wh2 * wv + wh1)

    # W update
    coeff_wu = (wh1 * wh2 + (nu + nonhydro) * wv)
    coeff_wv = nonhydro * (wh2 * wv - wh1)
    coeff_ww = (nu + nonhydro + wh2**2)

    Sol.rhou[i2] += -dt * thinv * coeff * denom * (coeff_uu * Dpx + coeff_uv * Dpy + coeff_uw * Dpz)
    Sol.rhov[i2] += -dt * thinv * coeff * denom * (coeff_vu * Dpx + coeff_vv * Dpy + coeff_vw * Dpz)
    if ndim == 3: Sol.rhow[i2] += -dt * thinv * coeff * denom * (coeff_wu * Dpx + coeff_wv * Dpy + coeff_ww * Dpz)

    # Sol.rhou[i2] += -dt * thinv * coeff * Dpx
    # Sol.rhov[i2] += -dt * thinv * coeff * Dpy
    # if ndim == 3: Sol.rhow[i2] += -dt * thinv * coeff * Dpz

    # Sol.rhou[...] = coeff_uu * 


    Sol.rhoX[i2] += - updt_chi * dt * dSdy * Sol.rhov[i2]


def operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt):
    g = ud.gravity_strength[1]
    Msq = ud.Msq
    Gammainv = th.Gammainv

    ndim = node.ndim
    nonhydro = ud.nonhydrostasy
    dy = elem.dy

    wh1, wv, wh2 = dt * ud.coriolis_strength

    ccenter = - ud.Msq * th.gm1inv / (dt**2)
    cexp = 2.0 - th.gamm

    igs = elem.igs

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
 
    strat = (mpv.HydroState_n.S0[y_idx1] - mpv.HydroState_n.S0[y_idx]) / dy

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

    nu = np.zeros_like(mpv.wcenter)
    nu[eindim] = -dt**2 * (g / Msq) * strat * Y

    setattr(mpv, 'nu_c', nu)
    nu = nu[eindim]

    denom = 1.0 / (wh1**2 + wh2**2 + (nu + nonhydro) * (wv**2 + 1))

    fimp = denom
    gimp = denom

    for dim in range(ndim):
        ## Assuming 2D vertical slice!
        if dim == 1:
            mpv.wplus[dim][eindim] = coeff #* gimp #* (wv**2 + 1.0)
        else:
            mpv.wplus[dim][eindim] = coeff #* fimp #* (wh1**2 + nu + nonhydro)

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
        is_wall = ud.bdry_type[dim] == BdryType.WALL or ud.bdry_type[dim] == BdryType.RAYLEIGH
        if is_wall:
            for direction in [-1,1]:
                wall_idx[dim] = (igs[dim]) * direction
                if direction == -1:
                    wall_idx[dim] -= 1
                wall_idx_tuple = tuple(wall_idx)
                rhs[wall_idx_tuple] *= factor


def grad_nodes(p2n, elem, node):
    ndim = node.ndim
    dx, dy, dz = node.dx, node.dy, node.dz

    kernels = []
    for dim in range(ndim):
        kernel = np.ones([2]*ndim)
        slc = [slice(None,)]*ndim
        slc[dim] = slice(0,1)
        kernel[tuple(slc)] *= -1.0
        kernels.append(kernel)

    dpdx = -0.5**(ndim-1) * signal.fftconvolve(p2n, kernels[0], mode='valid') / dx
    dpdy = -0.5**(ndim-1) * signal.fftconvolve(p2n, kernels[1], mode='valid') / dy if elem.iicy > 1 else 0.0
    dpdz = -0.5**(ndim-1) * signal.fftconvolve(p2n, kernels[2], mode='valid') / dz if (ndim == 3) else 0.0

    return dpdx, dpdy, dpdz



def divergence_nodes(rhs,elem,node,Sol,ud):
    ndim = elem.ndim
    igs = elem.igs
    dxyz = node.dxyz
    inner_idx = np.empty((ndim), dtype=object)

    for dim in range(ndim):
        is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
        inner_idx[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
    inner_idx_p1y = np.copy(inner_idx)
    inner_idx_p1y[1] = slice(1,-1)
        
    indices = [idx for idx in product([slice(0,-1),slice(1,None)], repeat = ndim)]
    signs = [sgn for sgn in product([1,-1], repeat = ndim)]
    inner_idx = tuple(inner_idx)
    inner_idx_p1y = tuple(inner_idx_p1y)


    if ud.bdry_type[1] == BdryType.WALL or ud.bdry_type[1] == BdryType.RAYLEIGH:
        Sol.rhou[:,:2,...] = 0.0 
        Sol.rhou[:,-2:,...] = 0.0

        Sol.rhov[:,:2,...] = 0.0 
        Sol.rhov[:,-2:,...] = 0.0 

        Sol.rhow[:,:2,...] = 0.0
        Sol.rhow[:,-2:,...] = 0.0

    Y = Sol.rhoY / Sol.rho

    Ux = np.diff(Sol.rhou * Y,axis=0) / elem.dx
    Ux = 0.5 * (Ux[:,:-1,...] + Ux[:,1:,...])

    Vy = np.diff(Sol.rhov * Y,axis=1) / elem.dy
    Vy = 0.5 * (Vy[:-1,...] + Vy[1:,...])

    if ndim == 3:
        Ux = -0.5 * (Ux[...,:-1] + Ux[...,1:])
        Vy = 0.5 * (Vy[...,:-1] + Vy[...,1:])

        Wz = np.diff(Sol.rhow * Y,axis=2) / elem.dz
        Wz = 0.5 * (Wz[:-1,...] + Wz[1:,...])
        Wz = 0.5 * (Wz[:,:-1,...] + Wz[:,1:,...])

        rhs[1:-1,1:-1,1:-1] = (Ux + Vy + Wz)
    else:
        rhs[1:-1,1:-1] = (Ux + Vy)

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
    rhs_n[inner_idx] = rhs[inner_idx] + 0.0 * rhs_hh
    return rhs_n
