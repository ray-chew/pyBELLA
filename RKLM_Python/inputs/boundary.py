"""
For more details on this module, refer to the write-up :ref:`boundary_handling`.
"""
from inputs.enum_bdry import BdryType

import numpy as np

from management.debug import find_nearest

def set_explicit_boundary_data(Sol, elem, ud, th, mpv, step=None):
    """
    In-place update of the ghost cells in :class:`management.variable.Vars` given the boundary conditions specified by :class:`inputs.user_data.UserDataInit`.

    Parameters
    ----------
    Sol : :class:`management.variable.Vars`
        Solution data container
    elem : :class:`discretization.kgrid.ElemSpaceDiscr`
        Cells grid
    ud : :class:`inputs.user_data.UserDataInit`
        Data container for the initial conditions
    th : :class:`physics.gas_dynamics.thermodynamic.ThemodynamicInit`
        Thermodynamic variables of the system
    mpv : :class:`physics.low_mach.mpv.MPV`
        Variables relating to the elliptic solver
    step : int, optional
        Current step

    """
    igs = elem.igs
    ndim = elem.ndim

    # if step parameter is not None, then we are in the advection directional Strang-splitting, where the array has already flipped, and we should only update the relevant boundaries, i.e. those in the direction of the current Strang-split-step.
    if step == None:
        dims = np.arange(ndim)
    else:
        dims = [1]

    for dim in dims:
        if step != None:
            current_step = step
        else:
            current_step = dim
        ghost_padding, idx = get_ghost_padding(ndim,dim,igs)

        if ud.gravity_strength[current_step] == 0.0:
            # Do this for the axes that do not have gravity.
            # Periodic BC.
            if ud.bdry_type[current_step] == BdryType.PERIODIC:
                set_boundary(Sol,ghost_padding,'wrap',idx,step=None)
            # Wall BC.
            else:
                set_boundary(Sol,ghost_padding,'symmetric',idx,step=None)

        else:
            # get current axis that has gravity.
            gravity_axis = dim
            
            direction = -1.
            offset = 0

            # get gravity strength specified in the user data file.
            g = ud.gravity_strength[gravity_axis]

            # for the number of ghost cells in the gravity axis...
            for side in ghost_padding[gravity_axis]:
                direction *= -1
                # print(direction)
                # loop through each of these ghost cells.
                for current_idx in np.arange(side)[::-1]:
                    nlast, nsource, nimage = get_gravity_padding(ndim,current_idx,direction,offset,elem)

                    # print(nlast, nsource, nimage)

                    Y_last = Sol.rhoY[nlast] / Sol.rho[nlast]
                    u = Sol.rhou[nsource] / Sol.rho[nsource]
                    w = Sol.rhow[nsource] / Sol.rho[nsource]
                    X = Sol.rhoX[nsource] / Sol.rho[nsource]

                    rhoYv_image = -Sol.rhov[nsource] * Sol.rhoY[nsource] / Sol.rho[nsource]

                    S = 1. / ud.stratification(elem.y[nimage[1]])

                    dpi = direction*(th.Gamma*g) * 0.5 * elem.dy * (1.0 / Y_last + S)

                    rhoY = ((Sol.rhoY[nlast]**th.gm1) + dpi)**th.gm1inv if ud.is_compressible == 1 else mpv.HydroState.rhoY0[nimage[1]]

                    rho = rhoY * S
                    p = rhoY**th.gamm

                    print(np.sign(direction))
                    if np.sign(direction) == 1:
                        v = rhoYv_image / rhoY
                    else:
                        v = -Sol.rhov[nsource] / Sol.rho[nsource]

                    Sol.rho[nimage] = rho

                    if np.sign(direction) == 1:
                        Sol.rhou[nimage] = Sol.rhou[nsource] * (Y_last * S)
                        Sol.rhow[nimage] = Sol.rhow[nsource] * (Y_last * S)
                    else:
                        Sol.rhou[nimage] = rho*u
                        Sol.rhov[nimage] = rho*v

                    Sol.rhow[nimage] = rho*w
                    # Sol.rhoe[nimage] = ud.rhoe(rho, u, v, w, p, ud, th)
                    Sol.rhoY[nimage] = rhoY
                    Sol.rhoX[nimage] = rho * X

                offset += 1


def set_boundary(Sol,pads,btype,idx,step=None):
    """
    Called by the function :func:`inputs.boundary.set_explicit_boundary_data`. Pads in-place the ghost cells for a given boundary type.

    Parameters
    ----------
    Sol : :class:`management.variable.Vars`
        Solution data container.
    pads : tuple
        A tuple containing the number of ghost cells to pad at each end.
    btype : string
        The type of boundary condition to pad. Currently supports:
            * `wrap` for periodic boundary conditions
            * `symmetric` for wall boundary conditions
            * `negative_symmetric` for wall boundary conditions, but with the signs flipped.
    idx : tuple
        A tuple containing the slice indices for the inner array, e.g. `(slice(2,-2),slice(2,-2))` for a 2D-array with 2 ghost cells for all edges.
    step : int, optional
        If we are in the advection routine with the flipped arrays according to the directional Strang-splitting, we want to pad the correct direction. `step=0`, pads the x-direction while `step=1` pads the y-direction.

    """
    Sol.rho[...] = np.pad(Sol.rho[idx],pads,btype)
    Sol.rhou[...] = np.pad(Sol.rhou[idx],pads,btype)
    if step == None or step == 0:
        Sol.rhou[...] = np.pad(Sol.rhou[idx],pads,btype)
        if btype == 'symmetric':
            Sol.rhov[...] = np.pad(Sol.rhov[idx],pads,negative_symmetric)
        else:
            Sol.rhov[...] = np.pad(Sol.rhov[idx],pads,btype)

    if step == 1:
        Sol.rhov[...] = np.pad(Sol.rhov[idx],pads,btype)
        if btype == 'symmetric':
            Sol.rhou[...] = np.pad(Sol.rhou[idx],pads,negative_symmetric)
        else:
            Sol.rhou[...] = np.pad(Sol.rhou[idx],pads,btype)

    Sol.rhow[...] = np.pad(Sol.rhow[idx],pads,btype)
    Sol.rhoe[...] = np.pad(Sol.rhoe[idx],pads,btype)
    Sol.rhoY[...] = np.pad(Sol.rhoY[idx],pads,btype)
    Sol.rhoX[...] = np.pad(Sol.rhoX[idx],pads,btype)


def negative_symmetric(vector,pad_width,iaxis,kwargs=None):
    """
    Taken from the reference:

    Parameters
    ----------
    vector : ndarray
        A rank 1 array already padded with zeros. Padded values are vector `[:iaxis_pad_width[0]] and vector[-iaxis_pad_width[1]:]`.
    iaxis_pad_width : tuple
        A 2-tuple of ints, `iaxis_pad_width[0]` represents the number of values padded at the beginning of vector where `iaxis_pad_width[1]` represents the number of values padded at the end of vector.
    iaxis : int
        The axis currently being calculated.
    kwargs : dict
        Any keyword arguments the function requires.

    References
    ----------
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.pad.html

    """
    if pad_width[1] > 0:
        sign = -1
        vector[:pad_width[0]] = sign * vector[pad_width[0]:2*pad_width[0]][::-1]
        vector[-pad_width[1]:] = sign * vector[-2*pad_width[1]:-pad_width[1]][::-1]
        return vector
    else: # axis must have length > 0 for padding
        return vector


def get_gravity_padding(ndim,cur_idx,direction,offset,elem,y_axs=None):
    """
    Parameters
    ----------
    ndim : int
        Number of dimensions.
    cur_idx : int
        The current index of the ghost cell in the gravity direction to be updated.
    direction : int
        Top of the domain, `direction=+1`, bottom of the domain, `direction=-1`.
    offset : int
        `offset=0`, index starts counting from 0,1.... `offset=1`, index starts counting from -1,-2,..., i.e. end-selection of the array.
    elem : :class:`discretization.kgrid.ElemSpaceDiscr`
        Cell grid.
    y_axs : int, optional
        `Default == None`. Specifies the direction of the gravity axis. If `None`, then direction is the the y-axis.

    """
    cur_i = np.copy(cur_idx)
    cur_idx += offset * ((elem.icy - 1) - 2*cur_idx)
    gravity_padding = [(slice(None))] * ndim
    if y_axs == None:
        # y_axs = ndim - 1
        y_axs = 1

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
    """
    For a given direction, return the number of ghost cells to pad the current direction, and the index slice of the inner array without the ghost cells.

    Parameters
    ----------
    ndim : int
        Number of dimensions for the problem.
    dim : int
        Current dimension to update
    igs : list
        A list of number of ghost cells in all dimensions, e.g. `[2,2,2]` for 2 ghost cells in the x, y, and z directions.

    Returns
    -------
    tuple
        Number of ghost cells in the current dimension at both edges.
    tuple
        Index slice of the inner domain of the array.

    """
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
    """
    Taken from the reference:

    Parameters
    ----------
    vector : ndarray
        A rank 1 array already padded with zeros. Padded values are vector `[:iaxis_pad_width[0]] and vector[-iaxis_pad_width[1]:]`.
    iaxis_pad_width : tuple
        A 2-tuple of ints, `iaxis_pad_width[0]` represents the number of values padded at the beginning of vector where `iaxis_pad_width[1]` represents the number of values padded at the end of vector.
    iaxis : int
        The axis currently being calculated.
    kwargs : dict
        Any keyword arguments the function requires.

    References
    ----------
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.pad.html

    """
    if all(pad_width) > 0:
        vector[:pad_width[0]+1], vector[-pad_width[1]-1:] = vector[-pad_width[1]-pad_width[1]-1:-pad_width[1]] , vector[pad_width[0]:pad_width[0]+pad_width[0]+1].copy()
    return vector


def check_flux_bcs(Lefts, Rights, elem, split_step, ud):
    igx = elem.igx
    igy = elem.igy

    if split_step == 1:
        if ud.bdry_type[split_step] == BdryType.WALL:

            left_inner = (slice(None),slice(igy,igy+1))
            left_ghost = (slice(None),slice(igy-1,igy))

            right_inner = (slice(None),slice(-igy-1,-igy))
            right_ghost = (slice(None),slice(-igy,-igy+1))

            rhou_wall = 0.
            # Lefts.rhou[left_ghost] = Rights.rhou[left_inner] = rhou_wall
            Lefts.rhoY[left_ghost] = Rights.rhoY[left_inner] = Rights.rho[left_inner] * ud.stratification(0.0)

            Lefts.rho[left_ghost] = Rights.rho[left_inner]
            Lefts.rhov[left_ghost] = Rights.rhov[left_inner]
            Lefts.rhow[left_ghost] = Rights.rhow[left_inner]
            Lefts.rhoe[left_ghost] = Rights.rhoe[left_inner]
            Lefts.rhoX[left_ghost] = Rights.rhoX[left_inner]
            # print(Lefts.rhoY[left_ghost])
            # print(Rights.rho[0,:])

            Rights.rho[right_ghost] = Lefts.rho[right_inner]
            Rights.rhou[right_ghost] = -1. * Lefts.rhou[right_inner]
            Rights.rhov[right_ghost] = Lefts.rhov[right_inner]
            Rights.rhow[right_ghost] = Lefts.rhow[right_inner]
            Rights.rhoe[right_ghost] = Lefts.rhoe[right_inner]
            Rights.rhoY[right_ghost] = Lefts.rhoY[right_inner]

            # print(Lefts.rhoY[:,-igx-2])
            
    else:
        if ud.bdry_type[split_step] == BdryType.WALL:

            assert(0) # INCOMPLETE!!!
            Lefts.rho[left_inner] = Rights.rho[:,igx-2]
            Lefts.rhou[left_inner] = -1. * Rights.rhou[:,igx-2]
            Lefts.rhov[left_inner] = Rights.rhov[:,igx-2]
            Lefts.rhow[left_inner] = Rights.rhow[:,igx-2]
            Lefts.rhoe[left_inner] = Rights.rhoe[:,igx-2]
            Lefts.rhoY[left_inner] = Rights.rhoY[:,igx-2]

            # print("#################### TRUE ########################")
            Rights.rho[right_ghost] = Lefts.rho[:,-igx-2]
            Rights.rhou[right_ghost] = -1. * Lefts.rhou[:,-igx-2]
            Rights.rhov[right_ghost] = Lefts.rhov[:,-igx-2]
            Rights.rhow[right_ghost] = Lefts.rhow[:,-igx-2]
            Rights.rhoe[right_ghost] = Lefts.rhoe[:,-igx-2]
            Rights.rhoY[right_ghost] = Lefts.rhoY[:,-igx-2]
            # print(Rights.rhoY[right_ghost])