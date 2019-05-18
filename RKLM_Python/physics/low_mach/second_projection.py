from inputs.enum_bdry import BdryType
from inputs.boundary import set_explicit_boundary_data
from scipy import signal
import numpy as np

def euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, dt, alpha_diff):
    nc = node.nc
    mpv.dp2_nodes = np.copy(mpv.p2_nodes)

    x_periodic = ud.bdry_type[0] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    z_periodic = ud.bdry_type[2] == BdryType.PERIODIC

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    operator_coefficients_nodes(elem, node, Sol, mpv, ud, th, dt)
    
    
    

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
    # igz = igs[2]
    nindim = np.empty((ndim),dtype='object')
    innerdim = np.copy(nindim)
    eindim = np.empty((ndim),dtype='object')

    for dim in range(ndim):
        is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
        nindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
        innerdim[dim] = slice(igs[dim],-igs[dim])
        eindim[dim] = slice(igs[dim]-is_periodic-1,-igs[dim]+is_periodic-2)

        if dim == 1:
            y_idx = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)
            right_idx = None if -igs[dim]+is_periodic == 0 else -igs[dim]+is_periodic
            y_idx1 = slice(igs[dim]-is_periodic+1, right_idx)

    strat = 2.0 * (mpv.HydroState_n.Y0[0,y_idx1] - mpv.HydroState_n.Y0[0,y_idx]) / (mpv.HydroState_n.Y0[0,y_idx1] + mpv.HydroState_n.Y0[0,y_idx]) / dy

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