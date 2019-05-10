from management.enumerator import HillShapes
from inputs.enum_bdry import BdryType

import numpy as np 

class InitializeBdry(object):
    def __init__(self, elem, ud):
        icx = elem.icx
        igx = elem.igx
        # igz = elem.igz
        self.wall_rhoYflux = np.zeros((icx-igx-igx))
        self.wall_slope = np.zeros((icx-igx-igx))
        self.wall_relative_slope = np.zeros((icx-igx-igx))

        if elem.ndim == 1:
            assert True, "This boundary condition setting makes no sense in 1D"

        elif elem.ndim == 2:
            self.wall_slope = slanted_wall_slope(elem.x[igx:-igx], ud)
            slope_sum = np.sum(np.abs(self.wall_slope))

            if slope_sum <= np.sqrt(np.finfo(np.float).eps):
                self.wall_relative_slope[:] = 0.0
            else:
                slope_sum_inv = 1.0 / slope_sum
                self.wall_relative_slope = slope_sum_inv * np.abs(self.wall_slope)

                slope_sum = np.sum(self.wall_relative_slope)

                print("relative slope sum = %e" %slope_sum)

            
def slanted_wall_slope_3d(x, size, ud):
    hill_height = ud.hill_height
    length_scale_inv = 1.0 / ud.hill_length_scale
    x_sc = x * length_scale_inv

    return np.repeat((-hill_height * 2.0 * x_sc / ((1.0 + x_sc * x_sc)**2) * length_scale_inv), size)
    
def slanted_wall_slope(x, ud):
    if (ud.hill_shape == HillShapes.SCHULTOW):
        kx = 2.0 * 2.0 * np.pi / (ud.xmax - ud.xmin)
        kz = np.sqrt(ud.Nsq / ud.wind_speed / ud.wind_speed + kx * kx)
        q = 0.25

        xi = kx * x
        y = 0.0

        for __ in range(10):
            y = q * np.cos(xi + y)

        yp = -q * np.sin(xi + y) / (1. + q * np.sin(xi + y))
        yp = kx * yp / kz
        return yp

    else:
        hill_height = ud.hill_height
        length_scale_inv = 1.0 / ud.hill_length_scale
        x_sc = x * length_scale_inv
        return (-hill_height * 2.0 * x_sc / ((1.0 + x_sc * x_sc)**2 * length_scale_inv)) * np.ones((x.shape))

def set_wall_rhoYflux(bdry, Sol0, mpv, elem, ud):
    igx = elem.igx
    # igy = elem.igy

    if elem.ndim == 1:
        print("wall flux in 1D makes no sense")
    elif elem.ndim == 2:
        is_x_periodic = True if (ud.bdry_type_min[0] == BdryType.PERIODIC) else False

        # idx_first_inner_j_row = (slice(igx,-igx),slice(0,igy))
        # print("idx_first_inner_j_row =", idx_first_inner_j_row)
        # print(Sol0.rho[idx_first_inner_j_row])
        bdry.wall_rhoYflux[...] = 0.0
        # bdry.wall_rhoYflux = wall_rhoYflux( \
        #                             elem.x[igx:-igx].reshape(-1,1), \
        #                             elem.y[igx:-igx].reshape(1,-1), \
        #                             Sol0.rhou[idx_first_inner_j_row]/ Sol0.rho[idx_first_inner_j_row], \
        #                             np.divide(Sol0.rhov[idx_first_inner_j_row],Sol0.rho[idx_first_inner_j_row]), \
        #                             mpv.HydroState_n.rhoY0[igy,igy], \
        #                             ud \
        #                             )
        # print(bdry.wall_rhoYflux)
        wall_flux_balance = np.sum(bdry.wall_rhoYflux)
        # correction for zero net flux
        bdry.wall_rhoYflux -= wall_flux_balance * bdry.wall_relative_slope

        # add the ghost cells back into wall_rhoYflux array: periodic or just zeroes?
        if (is_x_periodic):
            bdry.wall_rhoYflux = np.pad(bdry.wall_rhoYflux, 2, 'wrap')
        else:
            bdry.wall_rhoYflux = np.pad(bdry.wall_rhoYflux, 2, 'constant')

        flux_sum = np.sum(bdry.wall_rhoYflux[igx:-igx])
        print("wall flux sum = %e" %flux_sum)


def wall_rhoYflux(x,z,wind_speed_x,wind_speed_z,rhoY0,ud):
    return slanted_wall_slope(x,ud) * wind_speed_x * rhoY0

def set_explicit_boundary_data(Sol, elem, ud):
    igs = elem.igs
    ndim = elem.ndim

    for dim in range(ndim):
        ghost_padding, idx = get_ghost_padding(ndim,dim,igs)
        if ud.bdry_type[dim] == BdryType.PERIODIC:
            Sol.set_boundary(ghost_padding,'wrap',idx)
        else:
            Sol.set_boundary(ghost_padding,'symmetric',idx)


def get_ghost_padding(ndim,dim,igs):
    ghost_padding = [(0,0)] * ndim
    ghost_padding[dim] = (igs[dim],igs[dim])

    padded_idx = np.empty((ndim), dtype=object)
    for idim in range(ndim):
        padded_idx[idim] = slice(igs[idim],-igs[idim])
    padded_idx[dim] = slice(None)

    inner_domain = [slice(None)] * ndim
    inner_domain[dim] = slice(igs[dim],-igs[dim])

    return tuple(ghost_padding),  tuple(inner_domain)#tuple(padded_idx)

def bound(Sol, lambda_var, split_step, elem, ud):
    None

# def fancy(nsource,nlast,iimage,g,dh,th,bdry,Sol):
#     Y_last = Sol.rhoY[nlast] / Sol.rho[nlast]
#     v = Sol.rhov[nsource] / Sol.rho[nsource]
#     w = Sol.rhow[nsource] / Sol.rho[nsource]

#     Sol.rhoX = Sol.rhoX / Sol.rho

#     rhoYu_wall = bdry.wall_rhoYflux
#     rhoYu_image = 2.0 * rhoYu_wall - Sol.rhou[nsource] * Sol.rhoY[nsource] / Sol.rho[nsource]
#     S = 1. / ud.stratification(elem.x[iimage])

#     dpi = sign*(th.Gamma * g) * 0.5 * dh * (np.divide(1.0,Y_last) + S)
#     rhoY = (Sol.rhoY[nlast]**th.gm1 + dpi)**th.gm1inv
#     rho = rhoY * S
#     p = rhoY**th.gamm
#     if bd == 'top':
#         u = rhoYu_image / rhoY
#     else:
#         u = 2.0 * rhoYu_wall - Sol.rhou[nsource] / Sol.rho[nsource]

def periodic_plus_one(vector, pad_width, iaxis, kwargs=None):
    if all(pad_width) > 0:
        vector[:pad_width[0]+1], vector[-pad_width[1]-1:] = vector[-pad_width[1]-pad_width[1]-1:-pad_width[1]] , vector[pad_width[0]:pad_width[0]+pad_width[0]+1].copy()
    return vector

def set_ghostcells_p2(p,elem,ud):
    igs = elem.igs
    # idx = elem.inner_domain
    for dim in range(elem.ndim):
        ghost_padding, idx = get_ghost_padding(elem.ndim,dim,igs)
        if ud.bdry_type[dim] == BdryType.PERIODIC:
            p[...] = np.pad(p[idx],ghost_padding,'wrap')

        else: # WALL
            p[...] = np.pad(p[idx],ghost_padding,'symmetric')

def set_ghostnodes_p2(p,node,ud):
    igs = node.igs
    for dim in range(node.ndim):
        ghost_padding, idx = get_ghost_padding(node.ndim,dim,igs)

        if ud.bdry_type[dim] == BdryType.PERIODIC:
            p[...] = np.pad(p[idx], ghost_padding, periodic_plus_one)
        else: # ud.bdry_type[dim] == BdryType.WALL:
            p[...] = np.pad(p[idx], ghost_padding, 'symmetric')
