"""
For more details on this module, refer to the write-up :ref:`boundary_handling`.
"""

import numpy as np
from ....utils import options as opts
from .common import get_ghost_padding


def set_explicit_boundary_data(Sol, elem, ud, th, npf, step=None):
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
    th : :class:`physics.gas_dynamics.thermodynamic.init`
        Thermodynamic variables of the system
    npf : :class:`physics.low_mach.npf.MPV`
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
        dims = [ndim - 1]

    for dim in dims:
        if step is not None:
            current_step = step
        else:
            current_step = dim
        ghost_padding, idx = get_ghost_padding(ndim, dim, igs)

        if ud.gravity_strength[current_step] == 0.0:
            # Do this for the axes that do not have gravity.
            # Periodic BC.
            if ud.bdry_type[current_step] == opts.BdryType.PERIODIC:
                set_boundary(Sol, ghost_padding, "wrap", idx, step=None)
            # Wall BC.
            elif ud.bdry_type[current_step] == opts.BdryType.WALL:
                set_boundary(Sol, ghost_padding, "symmetric", idx, step=None)
            elif ud.bdry_type[current_step] == opts.BdryType.RAYLEIGH:
                assert 0, "Rayleigh boundary not defined on x-direction."

        else:
            # get current axis that has gravity.
            gravity_axis = dim

            direction = -1.0
            offset = 0

            # get gravity strength specified in the user data file.
            g = ud.gravity_strength[gravity_axis]

            # for the number of ghost cells in the gravity axis...
            for side in ghost_padding[gravity_axis]:
                direction *= -1
                # loop through each of these ghost cells.
                for current_idx in np.arange(side)[::-1]:
                    if step != None:
                        y_axs = ndim - 1
                    else:
                        y_axs = 1
                    nlast, nsource, nimage = get_gravity_padding(
                        ndim, current_idx, direction, offset, elem, y_axs=y_axs
                    )

                    Y_last = Sol.rhoY[nlast] / Sol.rho[nlast]

                    rhoYv_image = (
                        -Sol.rhov[nsource] * Sol.rhoY[nsource] / Sol.rho[nsource]
                    )

                    S = 1.0 / ud.stratification(elem.y[nimage[y_axs]])

                    if hasattr(ud, "ATMOSPHERIC_EXTENSION"):
                        dpi = (
                            npf.HydroState.p20[nimage[y_axs]]
                            - npf.HydroState.p20[nlast[y_axs]]
                        ) * ud.Msq
                    else:
                        dpi = (
                            direction
                            * (th.Gamma * g)
                            * 0.5
                            * elem.dy
                            * (1.0 / Y_last + S)
                        )

                    rhoY = (
                        ((Sol.rhoY[nlast] ** th.gm1) + dpi) ** th.gm1inv
                        if ud.is_compressible == 1
                        else npf.HydroState.rhoY0[nimage[y_axs]]
                    )

                    rho = rhoY * S

                    Y_source = Sol.rhoY[nsource] / Sol.rho[nsource]
                    Y_image = rhoY / rho

                    if hasattr(ud, "ATMOSPHERIC_EXTENSION"):
                        if direction > 0:  # if bottom boundary
                            v = Sol.rhov[nsource] * Y_source / Sol.rho[nsource] * rho
                        else:  # if top boundary
                            v = Sol.rhov[nsource] * Y_source

                        Th_slc = rhoY / rho / Y_last

                    else:
                        v = rhoYv_image / rhoY
                        Th_slc = 1.0

                    u = Sol.rhou[nsource] / Sol.rho[nsource]
                    w = Sol.rhow[nsource] / Sol.rho[nsource]
                    X = Sol.rhoX[nsource] / Sol.rho[nsource]

                    Sol.rho[nimage] = rho
                    Sol.rhou[nimage] = rho * u * Th_slc
                    if hasattr(ud, "ATMOSPHERIC_EXTENSION"):
                        Sol.rhov[nimage] = -v / Y_image
                    else:
                        Sol.rhov[nimage] = rho * v
                    Sol.rhow[nimage] = rho * w * Th_slc
                    Sol.rhoY[nimage] = rhoY
                    Sol.rhoX[nimage] = rho * X

                offset += 1


def set_boundary(Sol, pads, btype, idx, step=None):
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
    Sol.rho[...] = np.pad(Sol.rho[idx], pads, btype)

    if btype == "symmetric":
        Sol.rhov[...] = np.pad(Sol.rhov[idx], pads, negative_symmetric)
        Sol.rho[...] = np.pad(Sol.rho[idx], pads, "symmetric")
        Sol.rhou[...] = np.pad(Sol.rhou[idx], pads, "symmetric")
    elif btype == "constant":
        Sol.rho[...] = np.pad(Sol.rho[idx], pads, "symmetric")
        Sol.rhou[...] = np.pad(Sol.rhou[idx], pads, "symmetric")
        Sol.rhov[...] = np.pad(Sol.rhov[idx], pads, btype)
        Sol.rhow[...] = np.pad(Sol.rhow[idx], pads, "symmetric")
        btype = "symmetric"
    else:
        Sol.rhou[...] = np.pad(Sol.rhou[idx], pads, btype)
        Sol.rhov[...] = np.pad(Sol.rhov[idx], pads, btype)
        Sol.rhow[...] = np.pad(Sol.rhow[idx], pads, btype)

    Sol.rhoY[...] = np.pad(Sol.rhoY[idx], pads, btype)
    Sol.rhoX[...] = np.pad(Sol.rhoX[idx], pads, btype)


def negative_symmetric(vector, pad_width, iaxis, kwargs=None):
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
        vector[: pad_width[0]] = sign * vector[pad_width[0] : 2 * pad_width[0]][::-1]
        vector[-pad_width[1] :] = sign * vector[-2 * pad_width[1] : -pad_width[1]][::-1]
        return vector
    else:  # axis must have length > 0 for padding
        return vector


def get_gravity_padding(ndim, cur_idx, direction, offset, elem, y_axs=None):
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
    cur_idx += offset * ((elem.icy - 1) - 2 * cur_idx)
    gravity_padding = [slice(None)] * ndim
    if y_axs == None:
        # y_axs = ndim - 1
        y_axs = 1

    nlast = np.copy(gravity_padding)
    nlast[y_axs] = int(cur_idx + direction)

    nsource = np.copy(gravity_padding)
    nsource[y_axs] = int(
        offset * (elem.icy) + direction * (2 * elem.igy - (1 - offset) - cur_i)
    )

    nimage = np.copy(gravity_padding)
    nimage[y_axs] = int(cur_idx)
    return tuple(nlast), tuple(nsource), tuple(nimage)