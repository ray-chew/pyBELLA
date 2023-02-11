import numpy as np
from inputs.boundary import set_explicit_boundary_data, check_flux_bcs
from .recovery import recovery
from .numerical_flux import hll_solver

def advect(Sol, flux, dt, elem, odd, ud, th, mpv, node, label, writer = None):
    """
    Function that runs the advection routine with Strang-splitting. This function updates the `Sol` solution container with the advected solution in-place.

    Parameters
    ----------
    Sol : :py:class:`management.variable.Vars` 
        Solution data container.
    flux : :py:class:`management.variable.States`
        Fluxes data container
    dt : float
        Time-step size.
    elem : :py:class:`discretization.kgrid.ElemSpaceDiscr`
        Container for cell-grid properties.
    odd : int
        Is current step odd or even?
    ud : :py:class:`inputs.user_data.UserDataInit`
        Class container for initial conditions
    th : :py:class:`physics.gas_dynamics.thermodynamic.ThermodynamicInit`
        Class container for thermodynamic quantities.
    mpv : :py:class:`physics.low_mach.mpv.MPV`
        Container for Exner pressure.
    node : :py:class:`discretization.kgrid.NodeSpaceDiscr`
        Container for node-grid properties.
    label : string
        Tag label for the output array
    writer : :py:class:`management.io.io`, optional
        Writer class for I/O operations, by default None
    """
    # double strang sweep
    time_step = 0.5 * dt
    ndim = elem.ndim
    stage = 0

    # Sol.rhoX -= Sol.rho * mpv.HydroState.S0

    if (odd):
        for split in range(ndim):
            lmbda = time_step / elem.dxyz[split]
            Sol.flip_forward()
            if elem.isc[split] > 1: 
                explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)
    else:
        for i_split in range(ndim):
            split = elem.ndim - 1 - i_split
            lmbda = time_step / elem.dxyz[split]
            if elem.isc[split] > 1:
                explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv, [writer,node,label])
            Sol.flip_backward()

    stage = 1
    if (odd):
        for i_split in range(ndim):
            split = elem.ndim - 1 - i_split
            lmbda = time_step / elem.dxyz[split]
            if elem.isc[split] > 1:
                explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)
            Sol.flip_backward()
    else:
        for split in range(ndim):
            lmbda = time_step / elem.dxyz[split]
            Sol.flip_forward()
            if elem.isc[split] > 1:
                explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv)
            
    # Sol.rhoX += Sol.rho * mpv.HydroState.S0

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)



def explicit_step_and_flux(Sol, flux, lmbda, elem, split_step, stage, ud, th, mpv, writer=None, tag=None):
    """
    For each advection substep, solve the advection problem. For more details, see :ref:`advection_routine`. This function updates the solution `Sol` container in-place if a Strang-splitting is used, or returns the `flux` data container if a Runge-Kutta method is used.

    Parameters
    ----------
    Sol : :py:class:`management.variable.Vars` 
        Solution data container.
    flux : :py:class:`management.variable.States`
        Fluxes data container
    lmbda : float
        :math:`\\frac{dt}{dx}`, where :math:`dx` is the grid-size in the direction of the substep.
    elem : :py:class:`discretization.kgrid.ElemSpaceDiscr`
        Container for cell-grid properties.
    split_step : int
        Tracks the substep in the Strang-splitting.
    stage : int
        Tracks whether the substep order goes in x-y-z or z-y-x.
    ud : :py:class:`inputs.user_data.UserDataInit`
        Class container for initial conditions
    th : :py:class:`physics.gas_dynamics.thermodynamic.ThermodynamicInit`
        Class container for thermodynamic quantities.
    mpv : :py:class:`physics.low_mach.mpv.MPV`
        Container for Exner pressure.
    writer : :py:class:`management.io.io`, optional
        Writer class for I/O operations, by default None
    tag : str, optional
        If `rk` then the advection routine is solved by means of a first-order Runge-Kutta method, by default None

    Returns
    -------
    :py:class:`management.variable.States`
        `flux` data container.
    """
    set_explicit_boundary_data(Sol, elem, ud, th, mpv, step=split_step)

    Lefts, Rights = recovery(Sol, flux, lmbda, ud, th, elem, split_step, tag)

    # Lefts, Rights, u, Diffs, Ampls, Slopes = recovery(Sol, flux, lmbda, ud, th, elem, split_step, tag)

    # if writer is not None:
    #     writer[0].write_all(Sol,mpv,elem,writer[1],th,str(writer[2])+'_split_%i' %split_step)

    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'u',u)

    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Leftsu',Lefts.u)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Rightsu',Rights.u)

    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Leftsv',Lefts.v)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Rightsv',Rights.v)

    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Leftsw',Lefts.w)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Rightsw',Rights.w)

    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'LeftsrhoY',Lefts.rhoY)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'RightsrhoY',Rights.rhoY)

    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Diffsu',Diffs.u)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Diffsv',Diffs.v)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Diffsw',Diffs.w)

    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Amplsu',Ampls.u)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Amplsv',Ampls.v)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Amplsw',Ampls.w)

    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Slopesu',Slopes.u)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Slopesv',Slopes.v)
    # if writer is not None: writer[0].populate(str(writer[2])+'_split_%i' %split_step,'Slopesw',Slopes.w)

    # skipped check_flux_bcs for now; first debug other functions
    # check_flux_bcs(Lefts, Rights, elem, split_step, ud)

    flux = hll_solver(flux,Lefts,Rights,Sol, lmbda, ud, th)

    ndim = elem.ndim
    left_idx, right_idx = [slice(None)] * ndim, [slice(None)] * ndim
    right_idx[-1] = slice(1,None)
    left_idx[-1] = slice(0,-1)
    left_idx, right_idx = tuple(left_idx), tuple(right_idx)

    if tag != 'rk':
        Sol.rho += lmbda * (flux.rho[left_idx] - flux.rho[right_idx])
        Sol.rhou += lmbda * (flux.rhou[left_idx] - flux.rhou[right_idx])
        Sol.rhov += lmbda * (flux.rhov[left_idx] - flux.rhov[right_idx])
        Sol.rhow += lmbda * (flux.rhow[left_idx] - flux.rhow[right_idx])
        Sol.rhoe += lmbda * (flux.rhoe[left_idx] - flux.rhoe[right_idx])
        Sol.rhoX += lmbda * (flux.rhoX[left_idx] - flux.rhoX[right_idx])
        Sol.rhoY += lmbda * (flux.rhoY[left_idx] - flux.rhoY[right_idx])
    
    set_explicit_boundary_data(Sol, elem, ud, th, mpv, step=split_step)

    if tag == 'rk':
        return flux


def advect_rk(Sol, flux, dt, elem, odd, ud, th, mpv, node, label, writer = None):
    """
    Function that runs the advection routine with a first-order Runge-Kutta update. This function updates the `Sol` solution container with the advected solution in-place.

    Parameters
    ----------
    Sol : :py:class:`management.variable.Vars` 
        Solution data container.
    flux : :py:class:`management.variable.States`
        Fluxes data container
    dt : float
        Time-step size.
    elem : :py:class:`discretization.kgrid.ElemSpaceDiscr`
        Container for cell-grid properties.
    odd : int
        Is current step odd or even?
    ud : :py:class:`inputs.user_data.UserDataInit`
        Class container for initial conditions
    th : :py:class:`physics.gas_dynamics.thermodynamic.ThermodynamicInit`
        Class container for thermodynamic quantities.
    mpv : :py:class:`physics.low_mach.mpv.MPV`
        Container for Exner pressure.
    node : :py:class:`discretization.kgrid.NodeSpaceDiscr`
        Container for node-grid properties.
    label : string
        Tag label for the output array
    writer : :py:class:`management.io.io`, optional
        Writer class for I/O operations, by default None
        
    Attention
    ---------
    This function is not usually called unless commented out in the :py:meth:`management.data.time_update` routine.

    """
    # Do 1-stages Runge-Kutta.
    time_step = dt
    ndim = elem.ndim

    stage = 0
    # Get RK update
    for split in range(ndim):
        lmbda = time_step / elem.dxyz[split]
        Sol.flip_forward()
        if elem.isc[split] > 1: 
            flux[split] = explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv, tag='rk')

    ndim = elem.ndim
    left_idx, right_idx = [slice(None)] * ndim, [slice(None)] * ndim
    right_idx[-1] = slice(1,None)
    left_idx[-1] = slice(0,-1)
    left_idx, right_idx = tuple(left_idx), tuple(right_idx)

    for dim in range(ndim):
        lmbda = time_step / elem.dxyz[dim]
        Sol.flip_forward()
        Sol.rho += lmbda * (flux[dim].rho[left_idx] - flux[dim].rho[right_idx])
        Sol.rhou += lmbda * (flux[dim].rhou[left_idx] - flux[dim].rhou[right_idx])
        Sol.rhov += lmbda * (flux[dim].rhov[left_idx] - flux[dim].rhov[right_idx])
        Sol.rhow += lmbda * (flux[dim].rhow[left_idx] - flux[dim].rhow[right_idx])
        Sol.rhoe += lmbda * (flux[dim].rhoe[left_idx] - flux[dim].rhoe[right_idx])
        Sol.rhoX += lmbda * (flux[dim].rhoX[left_idx] - flux[dim].rhoX[right_idx])
        Sol.rhoY += lmbda * (flux[dim].rhoY[left_idx] - flux[dim].rhoY[right_idx])

        if dim == 1: # vertical axis
            updt = lmbda * (flux[dim].rhoX[left_idx] - flux[dim].rhoX[right_idx])
            setattr(Sol, 'pwchi', updt)

    set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    # stage = 1
    # for split in range(ndim):
    #     lmbda = dt / elem.dxyz[split]
    #     Sol.flip_forward()
    #     if elem.isc[split] > 1: 
    #         flux[split] = explicit_step_and_flux(Sol, flux[split], lmbda, elem, split, stage, ud, th, mpv, tag='rk')

    # Sol = deepcopy(Sol0)

    # for dim in range(ndim):
    #     lmbda = dt / elem.dxyz[dim]
    #     Sol.flip_forward()
    #     Sol.rho += lmbda * (flux[dim].rho[left_idx] - flux[dim].rho[right_idx])
    #     Sol.rhou += lmbda * (flux[dim].rhou[left_idx] - flux[dim].rhou[right_idx])
    #     Sol.rhov += lmbda * (flux[dim].rhov[left_idx] - flux[dim].rhov[right_idx])
    #     Sol.rhow += lmbda * (flux[dim].rhow[left_idx] - flux[dim].rhow[right_idx])
    #     Sol.rhoe += lmbda * (flux[dim].rhoe[left_idx] - flux[dim].rhoe[right_idx])
    #     Sol.rhoX += lmbda * (flux[dim].rhoX[left_idx] - flux[dim].rhoX[right_idx])
    #     Sol.rhoY += lmbda * (flux[dim].rhoY[left_idx] - flux[dim].rhoY[right_idx])
            
    # set_explicit_boundary_data(Sol, elem, ud, th, mpv)
