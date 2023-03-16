from inputs.enum_bdry import BdryType
from inputs.boundary import set_explicit_boundary_data, set_ghostnodes_p2, set_ghostcells_p2
from physics.low_mach.laplacian import stencil_9pt, stencil_27pt, stencil_hs, stencil_vs, stencil_9pt_numba_test, precon_diag_prepare
from scipy import signal
import numpy as np
from itertools import product
from numba import jit

from scipy.sparse.linalg import LinearOperator, bicgstab

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

    mpv.rhs[...] = divergence_nodes(mpv.rhs,elem,node,Sol,ud)
    if not hasattr(ud, 'LAMB_BDRY'): scale_wall_node_values(mpv.rhs, node, ud, 2.0)
    div = mpv.rhs

    if writer != None: writer.populate(str(label), 'rhs', div)

    rhoY = Sol.rhoY**(th.gamm - 2.0)
    dpidP_kernel = np.ones([2]*ndim)
    dpidP = (th.gm1 / ud.Msq) * signal.fftconvolve(rhoY, dpidP_kernel, mode='valid') / dpidP_kernel.sum()

    rhoYovG = Ginv * Sol.rhoY
    dbuoy = Sol.rhoY * (Sol.rhoX / Sol.rho)
    dpdx, dpdy, dpdz = grad_nodes(p2n, elem.ndim, node.dxyz)

    drhou = Sol.rhou - u0 * Sol.rho
    drhov = Sol.rhov - v0 * Sol.rho
    drhow = Sol.rhow - w0 * Sol.rho
    v = Sol.rhov / Sol.rho
    
    Sol.rhou = Sol.rhou - dt * ( rhoYovG * dpdx - corr_h2 * drhov + corr_v * drhow)

    Sol.rhov = Sol.rhov - dt * ( rhoYovG * dpdy + (g/Msq) * dbuoy * nonhydro - corr_h1 * drhow + corr_h2 * drhou) * (1 - ud.is_ArakawaKonor)

    Sol.rhow = Sol.rhow - dt * ( rhoYovG * dpdz - corr_v * drhou + corr_h1 * drhov) * (ndim == 3)

    Sol.rhoX = (Sol.rho * (Sol.rho / Sol.rhoY - S0c)) - dt * (v * dSdy) * Sol.rho

    dp2n[node.i1] -= dt * dpidP * div#[node.i1]

    weight = ud.compressibility
    mpv.p2_nodes += weight * dp2n

    set_ghostnodes_p2(mpv.p2_nodes,node, ud)
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)


def euler_backward_non_advective_expl_part(Sol, mpv, elem, dt, ud, th):
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[1]
    Msq = ud.Msq

    dbuoy = Sol.rhoY * (Sol.rhoX / Sol.rho)
    Sol.rhov = (nonhydro * Sol.rhov) - dt * (g/Msq) * dbuoy
    # Sol.rhov[np.where(Sol.rhov < 1e-15)] = 0.0

    Sol.mod_bg_wind(ud, -1.0)

    multiply_inverse_coriolis(Sol, Sol, mpv, ud, elem, elem, dt)

    Sol.mod_bg_wind(ud, +1.0)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)


total_iter = 0
total_calls = 0
def euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, dt, alpha_diff, Sol0 = None, writer = None, label=None):
    nc = node.sc
    rhs = np.zeros_like(mpv.p2_nodes)
    rhs = mpv.rhs

    # if elem.ndim == 2:
    #     # p2 = np.copy(mpv.p2_nodes[node.igx:-node.igx,node.igy:-node.igy])
    #     p2 = mpv.p2_nodes[node.p_isc]
    # elif elem.ndim == 3:
    #     p2 = np.copy(mpv.p2_nodes[1:-1,1:-1,1:-1])
    
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
        writer.populate(str(label),'wplusz',mpv.wplus[2]) if elem.ndim == 3 else writer.populate(str(label),'wplusz',np.zeros_like(mpv.wplus[0]))

    set_ghostnodes_p2(mpv.p2_nodes,node,ud)
    correction_nodes(Sol,elem,node,mpv,mpv.p2_nodes,dt,ud,th,0)
    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    rhs[...] = divergence_nodes(rhs,elem,node,Sol,ud)

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

    mpv.rhs[...] = rhs

    VS = True

    # prepare initial left-hand side and the laplacian stencil
    if elem.ndim == 2:
        Vec = mpv
        coriolis_params = multiply_inverse_coriolis(Vec, Sol, mpv, ud, elem, node, dt, attrs=('u', 'v', 'w'), get_coeffs = True)
        # lap = stencil_9pt(elem,node,mpv,Sol,ud,diag_inv,dt,coriolis_params)
        # sh = (ud.inx)*(ud.iny)

        diag_inv = precon_diag_prepare(mpv, elem, node, ud, coriolis_params)
        rhs *= diag_inv

        # diag_inv = np.ones_like(mpv.rhs)

        p2 = mpv.p2_nodes[node.i2].T
        lap = stencil_9pt_numba_test(mpv,node,coriolis_params,diag_inv, ud)
        sh = p2.shape[0] * p2.shape[1]

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
        # rhs_inner = rhs[node.igx:-node.igx,node.igy:-node.igy].ravel()
        # rhs_inner = rhs[1:-1,1:-1].ravel()
        rhs_inner = rhs[1:-1,1:-1].T.ravel()
        # rhs_inner = rhs.T.ravel()
    elif elem.ndim == 3 and elem.iicy > 1 and elem.iicz == 1:

        if not VS:
            rhs_inner = rhs[...,elem.igs[2]][node.igx:-node.igx,node.igy:-node.igy].ravel()
        if VS:
            rhs_inner = rhs[1:-1,1:-1,elem.igs[2]].ravel()

    elif elem.ndim == 3 and elem.icy - 2*elem.igs[1] > 2:
        rhs_inner = rhs[1:-1,1:-1,1:-1].ravel()
    else:
        rhs_inner = rhs[1:-1,elem.igs[1],1:-1].ravel()

    # p2, _ = bicgstab(lap,rhs_inner,tol=ud.tol,maxiter=ud.max_iterations,callback=counter)

    p2, _ = bicgstab(lap,rhs_inner,tol=ud.tol,maxiter=ud.max_iterations,callback=counter)
    # p2, _ = gmres(lap,rhs_inner,tol=ud.tol,maxiter=ud.max_iterations)
    # p2,info = bicgstab(lap,rhs.ravel(),x0=p2.ravel(),tol=1e-16,maxiter=6000,callback=counter)
    # print("Convergence info = %i, no. of iterations = %i" %(info,counter.niter))

    global total_calls, total_iter
    total_iter += counter.niter
    total_calls += 1
    print(counter.niter)
    print("Total calls to BiCGStab routine = %i, total iterations = %i" %(total_calls, total_iter))

    p2_full = np.zeros(nc).squeeze()
    if elem.ndim == 2:
        # p2_full[node.igx:-node.igx,node.igy:-node.igy] = p2.reshape(ud.inx,ud.iny)
        # p2_full[node.i2] = p2.reshape(rhs[node.i1].shape[0],rhs[node.i1].shape[1])
        p2_full[node.i2] = p2.reshape(rhs[node.i1].shape[1],rhs[node.i1].shape[0]).T
        # p2_full[node.i1] = p2.reshape(rhs.shape[1],rhs.shape[0]).T
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


def correction_nodes(Sol,elem,node,mpv,p,dt,ud,th,updt_chi):
    ndim = node.ndim
    Gammainv = th.Gammainv

    dSdy = mpv.HydroState_n.get_dSdy(elem, node)

    Dpx, Dpy, Dpz = grad_nodes(p, elem.ndim, node.dxyz)

    thinv = (Sol.rho / Sol.rhoY)

    Y = Sol.rhoY / Sol.rho
    coeff = (Gammainv * Sol.rhoY * Y)

    mpv.u = -dt * coeff * Dpx
    mpv.v = -dt * coeff * Dpy
    mpv.w = -dt * coeff * Dpz

    multiply_inverse_coriolis(mpv, Sol, mpv, ud, elem, elem, dt, attrs=['u', 'v', 'w'])

    Sol.rhou += thinv * mpv.u
    Sol.rhov += thinv * mpv.v
    Sol.rhow += thinv * mpv.w if ndim == 3 else 0.0
    Sol.rhoX += - updt_chi * dt * dSdy * Sol.rhov

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    assert True

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

    strat = mpv.HydroState_n.get_dSdy(elem, node)

    Y = Sol.rhoY / Sol.rho
    coeff = Gammainv * Sol.rhoY * Y

    nu = np.zeros_like(mpv.wcenter)
    nu = -dt**2 * (g / Msq) * strat * Y

    setattr(mpv, 'nu_c', nu)
    nu = nu

    # coeff[:,0] = coeff[:,2]
    # coeff[:,1] = coeff[:,2]
    # coeff[:,-1] = coeff[:,-3]
    # coeff[:,-2] = coeff[:,-3]
    # coeff[:,:2] = coeff[:,-2:] = 0.0
    # coeff[0,:] = coeff[-1,:] = 0.0

    for dim in range(ndim):
        ## Assuming 2D vertical slice!
        if dim == 1:
            mpv.wplus[dim][...] = coeff #* gimp #* (wv**2 + 1.0)
        else:
            mpv.wplus[dim][...] = coeff #* fimp #* (wh1**2 + nu + nonhydro)

    kernel = np.ones([2] * ndim)

    mpv.wcenter = (ccenter * signal.fftconvolve(Sol.rhoY**cexp,kernel,mode='valid') / kernel.sum())

    tmp_wplus = signal.fftconvolve(mpv.wplus[0],kernel,mode='valid') / kernel.sum()

    # set_ghostcells_p2(mpv.wplus[0], elem, ud)
    # set_ghostcells_p2(mpv.wplus[1], elem, ud)
    # set_ghostnodes_p2(mpv.wcenter, node, ud, igs=(1,1))
    # mpv.wcenter[:,0] = mpv.wcenter[:,1]
    # mpv.wcenter[:,-1] = mpv.wcenter[:,-2]

    assert True
    if not hasattr(ud, 'LAMB_BDRY'): scale_wall_node_values(mpv.wcenter, node, ud)


# def operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt):
#     g = ud.gravity_strength[1]
#     Msq = ud.Msq
#     Gammainv = th.Gammainv

#     ndim = node.ndim
#     nonhydro = ud.nonhydrostasy
#     dy = elem.dy

#     wh1, wv, wh2 = dt * ud.coriolis_strength

#     ccenter = - ud.Msq * th.gm1inv / (dt**2)
#     cexp = 2.0 - th.gamm

#     igs = elem.igs

#     nindim = np.empty((ndim),dtype='object')
#     innerdim = np.empty((ndim),dtype='object')
#     innerdim1 = np.empty((ndim),dtype='object')
#     eindim = np.empty((ndim),dtype='object')

#     for dim in range(ndim):
#         is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
#         nindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
#         innerdim[dim] = slice(igs[dim],-igs[dim])
#         eindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)

#         if dim == 1:
#             y_idx = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)
#             right_idx = None if -igs[dim]+is_periodic == 0 else -igs[dim]+is_periodic
#             y_idx1 = slice(igs[dim]-is_periodic+1, right_idx)

#         innerdim1[dim] = slice(igs[dim]-1, (-igs[dim]+1))
 
#     strat = (mpv.HydroState_n.S0[y_idx1] - mpv.HydroState_n.S0[y_idx]) / dy

#     nindim = tuple(nindim)
#     eindim = tuple(eindim)
#     innerdim = tuple(innerdim)
#     innerdim1 = tuple(innerdim1)

#     for dim in range(0,elem.ndim,2):
#         is_periodic = ud.bdry_type[dim] != BdryType.PERIODIC
#         strat = np.expand_dims(strat, dim)
#         strat = np.repeat(strat, elem.sc[dim]-int(2*is_periodic+igs[dim]), axis=dim)

#     Y = Sol.rhoY[nindim] / Sol.rho[nindim]
#     coeff = Gammainv * Sol.rhoY[nindim] * Y

#     nu = np.zeros_like(mpv.wcenter)
#     nu[eindim] = -dt**2 * (g / Msq) * strat * Y

#     setattr(mpv, 'nu_c', nu)
#     nu = nu[eindim]

#     denom = 1.0 / (wh1**2 + wh2**2 + (nu + nonhydro) * (wv**2 + 1))

#     fimp = denom
#     gimp = denom

#     for dim in range(ndim):
#         ## Assuming 2D vertical slice!
#         if dim == 1:
#             mpv.wplus[dim][eindim] = coeff #* gimp #* (wv**2 + 1.0)
#         else:
#             mpv.wplus[dim][eindim] = coeff #* fimp #* (wh1**2 + nu + nonhydro)

#     kernel = np.ones([2] * ndim)

#     mpv.wcenter[innerdim] = ccenter * signal.fftconvolve(Sol.rhoY[innerdim1]**cexp,kernel,mode='valid') / kernel.sum()

#     scale_wall_node_values(mpv.wcenter, node, ud)


def scale_wall_node_values(rhs, node, ud, factor=.5):
    # if factor < 1.0:
    # if factor < 1.0:
    #     factor = 1.0
    #     rhs[:,1] *= factor
    #     rhs[:,-2] *= factor

    #     # rhs[:,0] *=  1.0# factor
    #     # rhs[:,-1] *= 1.0# factor
    # rhs[:,0] = rhs[:,2] * factor
    # rhs[:,-1] = rhs[:,-3] * factor
    # else:
    #     factor = 1.0
    #     rhs[:,:2] *= factor
    #     rhs[:,-2:] *= factor

    ndim = node.ndim
    igs = node.igs

    wall_idx = np.empty((ndim), dtype=object)
    for dim in range(ndim):
        wall_idx[dim] = slice(igs[dim],-igs[dim])

    for dim in range(ndim):
        is_wall = ud.bdry_type[dim] == BdryType.WALL or ud.bdry_type[dim] == BdryType.RAYLEIGH
        if is_wall:
            for direction in [-1,1]:
                wall_idx[dim] = (igs[dim]-1) * direction
                if direction == -1:
                    wall_idx[dim] -= 1
                wall_idx_tuple = tuple(wall_idx)
                rhs[wall_idx_tuple] *= factor

        # is_rayleigh = ud.bdry_type[dim] == BdryType.RAYLEIGH
        # if is_rayleigh:
        #     for direction in [1]:
        #         wall_idx[dim] = (igs[dim]-1) * direction
        #         # if direction == -1:
        #         #     wall_idx[dim] -= 1
        #         wall_idx_tuple = tuple(wall_idx)
        #         rhs[wall_idx_tuple] *= factor



def grad_nodes_fft(p2n, elem, node):
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

# @jit(nopython=True, nogil=False, cache=True)
def grad_nodes(p, ndim, dxy):
    dx, dy, dz = dxy

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

    Dpx *= 0.5**(ndim-1) / dx
    Dpy *= 0.5**(ndim-1) / dy
    Dpz *= 0.5**(ndim-1) / dz

    return Dpx, Dpy, Dpz



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

    if not hasattr(ud, 'LAMB_BDRY'):
        if ud.bdry_type[1] == BdryType.WALL or ud.bdry_type[1] == BdryType.RAYLEIGH:
            Sol.rhou[:,:2,...] = 0.0 
            Sol.rhov[:,:2,...] = 0.0  
            Sol.rhow[:,:2,...] = 0.0

        # if ud.bdry_type[1] == BdryType.WALL:
            Sol.rhou[:,-2:,...] = 0.0
            Sol.rhov[:,-2:,...] = 0.0
            Sol.rhow[:,-2:,...] = 0.0

    # if ud.bdry_type[1] == BdryType.RAYLEIGH:
    #     u_last = Sol.rhou[:,-3,...] / Sol.rho[:,-3,...]
    #     v_last = Sol.rhov[:,-3,...] / Sol.rho[:,-3,...]

    #     Sol.rhou[:,-2,...] = u_last * Sol.rho[:,-2,...]
    #     Sol.rhov[:,-2,...] = -v_last * Sol.rho[:,-2,...]

    #     Sol.rhou[:,-1,...] = u_last * Sol.rho[:,-1,...]
    #     Sol.rhov[:,-1,...] = -(Sol.rhov[:,-4,...] / Sol.rho[:,-4,...]) * Sol.rho[:,-4,...]

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
        rhs = (Ux + Vy)

    rhs_max = np.max(rhs[inner_idx]) if np.max(rhs[inner_idx]) > 0 else 0
    return rhs


def rhs_from_p_old(rhs,node,mpv):
    igs = node.igs
    ndim = node.ndim

    assert ndim != 1, "Not implemented for 1D"
    
    # inner_idx = np.empty((ndim), dtype=object)
    # for dim in range(ndim):
    #     inner_idx[dim] = slice(igs[dim],-igs[dim])
    
    # inner_idx = tuple(inner_idx)
    # rhs_n = np.zeros_like(rhs)
    # rhs_hh = mpv.wcenter[inner_idx] * mpv.p2_nodes[inner_idx]
    # rhs_n[inner_idx] = rhs[inner_idx] + 0.0 * rhs_hh

    inner_idx = np.empty((ndim), dtype=object)
    for dim in range(ndim):
        inner_idx[dim] = slice(igs[dim],-igs[dim])
    
    inner_idx = tuple(inner_idx)
    rhs_n = np.zeros_like(rhs)
    rhs_hh = mpv.wcenter * mpv.p2_nodes[node.i1]
    rhs_n = rhs + 0.0 * rhs_hh
    return rhs_n


def multiply_inverse_coriolis(Vec, Sol, mpv, ud, elem, node, dt, attrs=('rhou', 'rhov', 'rhow'), get_coeffs = False):
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[1]
    Msq = ud.Msq

    wh1, wv, wh2 = dt * ud.coriolis_strength

    strat = mpv.HydroState_n.get_dSdy(elem, node)

    Y = Sol.rhoY / Sol.rho
    nu = -dt**2 * (g / Msq) * strat * Y

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

    VecU = getattr(Vec, attrs[0])
    VecV = getattr(Vec, attrs[1])
    VecW = getattr(Vec, attrs[2])

    # Do the updates
    U = denom * (coeff_uu * VecU + coeff_uv * VecV + coeff_uw * VecW)
    V = denom * (coeff_vu * VecU + coeff_vv * VecV + coeff_vw * VecW)
    W = denom * (coeff_wu * VecU + coeff_wv * VecV + coeff_ww * VecW)

    VecU[...] = U
    VecV[...] = V
    VecW[...] = W

    i1 = node.i1
    if get_coeffs:
        # coriolis_parameters = ((coeff_uu * denom)[i1].reshape(-1,), (coeff_vv * denom)[i1].reshape(-1,), (coeff_uv * denom)[i1].reshape(-1,), (coeff_vu * denom)[i1].reshape(-1,))
        coriolis_parameters = ((coeff_uu * denom).T, (coeff_vv * denom).T, (coeff_uv * denom).T, (coeff_vu * denom).T)
        return coriolis_parameters