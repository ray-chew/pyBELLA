from management.enumerator import HillShapes
from inputs.enum_bdry import BdryType

import numpy as np

def set_explicit_boundary_data(Sol, elem, ud, th, mpv):
    igs = elem.igs
    ndim = elem.ndim

    for dim in range(ndim):
        ghost_padding, idx = get_ghost_padding(ndim,dim,igs)

        if ud.gravity_strength[dim] == 0.0:
            if ud.bdry_type[dim] == BdryType.PERIODIC:
                set_boundary(Sol,ghost_padding,'wrap',idx)
            else:
                set_boundary(Sol,ghost_padding,'symmetric',idx)

        else:
            # recursive updating of array - for loop cannot be avoided....?
            # assumption: gravity always acts in the y-axis
            gravity_axis = 1
            direction = -1.
            offset = 0

            g = ud.gravity_strength[gravity_axis]

            for side in ghost_padding[gravity_axis]:
                direction *= -1
                for current_idx in np.arange(side)[::-1]:

                    nlast, nsource, nimage = get_gravity_padding(ndim,current_idx,direction,offset,elem)

                    Y_last = Sol.rhoY[nlast] / Sol.rho[nlast]
                    u = Sol.rhou[nsource] / Sol.rho[nsource]
                    w = Sol.rhow[nsource] / Sol.rho[nsource]
                    X = Sol.rhoX[nsource] / Sol.rho[nsource]

                    rhoYv_image = Sol.rhov[nsource] * Sol.rhoY[nsource] / Sol.rho[nsource]

                    S = 1. / ud.stratification(elem.y[nimage[1]])

                    dpi = direction*(th.Gamma*g) * 0.5 * elem.dy * (1.0 / Y_last + S)
                    rhoY = ((Sol.rhoY[nlast]**th.gm1) + dpi)**th.gm1inv if ud.is_compressible == 1 else mpv.HydroState.rhoY0[nimage]

                    rho = rhoY * S
                    p = rhoY**th.gamm
                    if np.sign(direction) == 1:
                        v = rhoYv_image / rhoY
                    else:
                        v = -Sol.rhov[nsource] / Sol.rho[nsource]

                    Sol.rho[nimage] = rho
                    Sol.rhou[nimage] = rho*u
                    Sol.rhov[nimage] = rho*v
                    Sol.rhow[nimage] = rho*w
                    Sol.rhoe[nimage] = ud.rhoe(rho, u, v, w, p, ud, th)
                    Sol.rhoY[nimage] = rhoY
                    Sol.rhoX[nimage] = rho * X

                offset += 1


def set_boundary(Sol,pads,btype,idx):
    Sol.rho[...] = np.pad(Sol.rho[idx],pads,btype)
    if btype == 'symmetric':
        Sol.rhou[...] = np.pad(Sol.rhou[idx],pads,negative_symmetric)
    else:
        Sol.rhou[...] = np.pad(Sol.rhou[idx],pads,btype)
    Sol.rhou[...] = np.pad(Sol.rhou[idx],pads,btype)
    Sol.rhov[...] = np.pad(Sol.rhov[idx],pads,btype)
    Sol.rhow[...] = np.pad(Sol.rhow[idx],pads,btype)
    Sol.rhoe[...] = np.pad(Sol.rhoe[idx],pads,btype)
    Sol.rhoY[...] = np.pad(Sol.rhoY[idx],pads,btype)
    Sol.rhoX[...] = np.pad(Sol.rhoX[idx],pads,btype)


def negative_symmetric(vector,pad_width,iaxis,kwargs=None):
    if pad_width[1] > 0:
        sign = 1
        vector[:pad_width[0]] = sign * vector[pad_width[0]:2*pad_width[0]][::-1]
        vector[-pad_width[1]:] = sign * vector[-2*pad_width[1]:-pad_width[1]][::-1]
        return vector
    else: # axis must have length > 0 for padding
        return vector


def get_gravity_padding(ndim,cur_idx,direction,offset,elem):
    cur_i = np.copy(cur_idx)
    cur_idx += offset * ((elem.icy - 1) - 2*cur_idx)
    gravity_padding = [(slice(None))] * ndim
    y_axs = ndim - 1

    nlast = np.copy(gravity_padding)
    nlast[y_axs] = int(cur_idx + direction)

    nsource = np.copy(gravity_padding)
    nsource[y_axs] = int(offset*(elem.icy) + direction * (2 * elem.igy - (1 - offset) - cur_i))

    nimage = np.copy(gravity_padding)
    nimage[y_axs] = int(cur_idx)
    return tuple(nlast), tuple(nsource), tuple(nimage)


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
            p[...] = np.pad(p[idx], ghost_padding, 'reflect')


def get_ghost_padding(ndim,dim,igs):
    ghost_padding = [(0,0)] * ndim
    ghost_padding[dim] = (igs[dim],igs[dim])

    padded_idx = np.empty((ndim), dtype=object)
    for idim in range(ndim):
        padded_idx[idim] = slice(igs[idim],-igs[idim])
    padded_idx[dim] = slice(None)

    inner_domain = [slice(None)] * ndim
    inner_domain[dim] = slice(igs[dim],-igs[dim])

    return tuple(ghost_padding),  tuple(inner_domain)


def periodic_plus_one(vector, pad_width, iaxis, kwargs=None):
    if all(pad_width) > 0:
        vector[:pad_width[0]+1], vector[-pad_width[1]-1:] = vector[-pad_width[1]-pad_width[1]-1:-pad_width[1]] , vector[pad_width[0]:pad_width[0]+pad_width[0]+1].copy()
    return vector