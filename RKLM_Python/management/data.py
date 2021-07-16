from discretization.kgrid import Grid, ElemSpaceDiscr, NodeSpaceDiscr

# dependencies of the atmospheric flow solver
import inputs.boundary as boundary
from management.variable import States, Vars
from discretization import kgrid
from physics.gas_dynamics.thermodynamic import ThermodynamicInit
from physics.gas_dynamics.numerical_flux import recompute_advective_fluxes
from physics.gas_dynamics.explicit import advect, advect_rk
from physics.gas_dynamics.eos import nonhydrostasy, compressibility, synchronise_variables, is_compressible, is_nonhydrostatic
from physics.gas_dynamics.gas_dynamics import dynamic_timestep
from physics.low_mach.second_projection import euler_backward_non_advective_impl_part, euler_backward_non_advective_expl_part, euler_forward_non_advective
from inputs.enum_bdry import BdryType
from physics.low_mach.mpv import MPV

# blending module
from data_assimilation import blending

import numpy as np
from copy import deepcopy
from scipy import signal
import h5py

from termcolor import colored

def data_init(ud):
    """
    Helper function to initialise the `elem` and `node` grids, corresponding to the cell and node grids, from a given user iniital data file.

    Parameters
    ----------
    ud : :class:`inputs.user_data.UserDataInit`
        Data container for the initial conditions.

    Returns
    -------
    elem : :class:`discretization.kgrid.ElemSpaceDiscr`
        Cells grid.
    node : :class:`discretization.kgrid.NodeSpaceDiscr`
        Nodes grid.

    """
    inx = ud.inx
    iny = ud.iny
    inz = ud.inz
    x0 = ud.xmin
    x1 = ud.xmax
    y0 = ud.ymin
    y1 = ud.ymax
    z0 = ud.zmin
    z1 = ud.zmax

    grid = Grid(inx,iny,inz,x0,x1,y0,y1,z0,z1)

    elem = ElemSpaceDiscr(grid)
    node = NodeSpaceDiscr(grid)

    return elem, node

def time_update(Sol,flux,mpv,t,tout,ud,elem,node,steps,th,bld=None,writer=None,debug=False):
    """
    For more details, refer to the write-up :ref:`time-stepping`.

    Does a time-step for the atmospheric solver.

    Parameters
    ----------
    Sol : :class:`management.variable.Vars`
        Solution data container.
    flux : :class:`management.variable.States`
        Data container for the fluxes.
    mpv : :class:`physics.low_mach.mpv.MPV`
        Variables relating to the elliptic solver.
    t : float
        Current time
    tout : float
        Next output time
    ud : :class:`inputs.user_data.UserDataInit`
        Data container for the initial conditions
    elem : :class:`discretization.kgrid.ElemSpaceDiscr`
        Cells grid.
    node : :class:`discretization.kgrid.NodeSpaceDiscr`
        Nodes grid.
    step : int
        Current step.
    th : :class:`physics.gas_dynamics.thermodynamic.ThermodynamicInit`
        Thermodynamic variables of the system
    bld : :class:`data_assimilation.blending.Blend()`
        Blending class used to initalise interface blending methods.
    writer : :class:`management.io.io`, optional
        `default == None`. If given, output after each time-step will be written in the hdf5 format.  
    debug : boolean, optional
        `default == False`. If `True`, then writer will output `Sol`:
            1. before flux calculation 
            2. before advection routine 
            3. after advection routine 
            4. after explicit solver
            5. after implicit solver
        
        during both the half-step for the prediction of advective flux and the full-step.

    Returns
    -------
    list
        A list of `[Sol,flux,mpv,[window_step,step]]` data containers at time `tout`.
    """

    window_step = steps[0]
    step = steps[1]
    swe_to_lake = False
    if 'best' not in ud.aux:
        test_hydrob = True
    else:
        test_hydrob = False

    while ((t < tout) and (step < ud.stepmax)):
        boundary.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

        label = '%.3d' %step

        if step == 0 and writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_ic')

        dt, cfl, cfl_ac = dynamic_timestep(Sol,t,tout,elem,ud,th, step)

        if 'CFLfixed' in ud.aux:
            if step < 2 : dt = 21.69 / ud.t_ref


        c1 = step == 0 and ud.is_nonhydrostatic == 0 and bld is not None and 'imbal' in ud.aux
        c2 = step == 0 and ud.is_nonhydrostatic == 1 and ud.initial_blending == True and bld is not None and 'imbal' in ud.aux

        if c1 or c2:
            print(colored("nonhydrostatic to hydrostatic conversion...", 'blue'))
            ud.is_nonhydrostatic = 0
            if test_hydrob == False:
                dt *= 0.5

        ######################################################
        # Blending : Do blending before timestep
        ######################################################
        swe_to_lake, Sol, mpv, t = blending.blending_before_timestep(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt, swe_to_lake, debug)

        ud.is_nonhydrostatic = is_nonhydrostatic(ud,window_step)
        ud.nonhydrostasy = nonhydrostasy(ud,t,window_step)

        if ud.continuous_blending == True or ud.initial_blending == True:
            if window_step >= 0:
                print("step = %i, window_step = %i" %(step,window_step))
            else:
                print("step = %i, window_step = %f" %(step,window_step))
        print("is_compressible = %i, is_nonhydrostatic = %i" %(ud.is_compressible, ud.is_nonhydrostatic))
        print("compressibility = %.3f, nonhydrostasy = %.3f" %(ud.compressibility,ud.nonhydrostasy))
        print("-------")

        Sol0 = deepcopy(Sol)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_flux')
        
        recompute_advective_fluxes(flux, Sol)

        if debug == True: writer.populate(str(label)+'_before_advect','rhoYu',flux[0].rhoY)
        if debug == True: writer.populate(str(label)+'_before_advect','rhoYv',flux[1].rhoY)
        if debug == True and elem.ndim == 3: writer.populate(str(label)+'_before_advect','rhoYw',flux[2].rhoY)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_advect')

        advect(Sol, flux, 0.5*dt, elem, step%2, ud, th, mpv, node, str(label)+'_half', writer)
        # advect_rk(Sol, flux, 0.5*dt, elem, step%2, ud, th, mpv, node, str(label)+'_half', writer)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_advect')

        if writer is not None: writer.populate(str(label)+'_after_full_step','p2_nodes0',mpv.p2_nodes)

        mpv.p2_nodes0[...] = mpv.p2_nodes

        euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaexp')
        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 1.0, label=str(label+'_after_ebnaimp'), writer=writer)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaimp')

        # if test_hydrob == True and writer is not None and step==0:
        #     writer.write_all(Sol,mpv,elem,node,th,str(label)+'_half')

        flux_half_new = deepcopy(flux)

        recompute_advective_fluxes(flux, Sol)

        if debug == True: writer.populate(str(label)+'_after_half_step','rhoYu',flux[0].rhoY)
        if debug == True: writer.populate(str(label)+'_after_half_step','rhoYv',flux[1].rhoY)
        if debug == True and elem.ndim == 3: writer.populate(str(label)+'_after_half_step','rhoYw',flux[2].rhoY)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_half_step')

        Sol_half_new = deepcopy(Sol)
        mpv_half_new = deepcopy(mpv)
        rho_half = np.copy(Sol.rho)
        rhou_half = np.copy(Sol.rhou)
        rhov_half = np.copy(Sol.rhov)
        rhow_half = np.copy(Sol.rhow)
        rhoX_half = np.copy(Sol.rhoX)
        rhoY_half = np.copy(Sol.rhoY)
        p2_nodes_half = np.copy(mpv.p2_nodes)

        if test_hydrob == True and writer is not None and step==0:
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_half')
        #     print(dt*0.5)

        # takes care of the non-hydrostatic, compressible case
        if ud.is_compressible == 1 and ud.is_nonhydrostatic == 1:
            mpv.p2_nodes[...] = mpv.p2_nodes0
            Sol = deepcopy(Sol0)
        # takes care of the hydrostatic case
        elif ud.is_nonhydrostatic == 0:
            # if step == 1 :
                # Sol.rho = (Sol.rho_half + Sol.rho) * 0.5
                # Sol.rhou = (Sol.rhou_half + Sol.rhou) * 0.5
                # Sol.rhov = (Sol.rhov_half + Sol.rhov) * 0.5
                # Sol.rhow = (Sol.rhow_half + Sol.rhow) * 0.5
                # Sol.rhoX = (Sol.rhoX_half + Sol.rhoX) * 0.5
                # Sol.rhoY = (Sol.rhoY_half + Sol.rhoY) * 0.5
                # v = (Sol.rhov_half / Sol.rho_half + Sol.rhov / Sol.rho) * 0.5
                # mpv.p2_nodes = (mpv.p2_nodes_half + mpv.p2_nodes) * 0.5
                # v = rhov_half / rho_half
                # Sol = deepcopy(Sol0)
                # Sol.rhov = (Sol.rhov_half + rhov_half) * 0.5
                # Sol.rhov = Sol.rho * v
                # mpv.p2_nodes[...] = mpv.p2_nodes0
                # None
                # Sol = deepcopy(Sol_tu)
                # mpv = deepcopy(mpv_tu)
                # None
            # else:
            Sol = deepcopy(Sol0)
            mpv.p2_nodes[...] = mpv.p2_nodes0
            
        # takes care of the pseudo-incompressible case
        elif ud.is_compressible == 0:
            Sol = deepcopy(Sol0)

        Sol.rhov0 = np.copy(Sol.rhov)
        Sol.rho_half = rho_half
        Sol.rhou_half = rhou_half
        Sol.rhov_half = rhov_half
        Sol.rhow_half = rhow_half
        Sol.rhoX_half = rhoX_half
        Sol.rhoY_half = rhoY_half
        mpv.p2_nodes_half = p2_nodes_half

        # Sol.rho_half = Sol.rho
        # Sol.rhou_half = Sol.rhou
        # Sol.rhov_half = Sol.rhov
        # Sol.rhow_half = Sol.rhow
        # Sol.rhoX_half = Sol.rhoX
        # Sol.rhoY_half = Sol.rhoY
        # mpv.p2_nodes_half = mpv.p2_nodes
        Sol_half_old = deepcopy(Sol_half_new)
        flux_half_old = deepcopy(flux_half_new)
        mpv_half_old = deepcopy(mpv_half_new)

        # mpv.p2_nodes[...] = ud.compressibility * mpv.p2_nodes0 + (1.0-ud.compressibility) * mpv.p2_nodes

        euler_forward_non_advective(Sol, mpv, elem, node, 0.5*dt, ud, th, writer=writer, label=str(label)+'_after_efna')

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_efna')
 
        advect(Sol, flux, dt, elem, step%2, ud, th, mpv, node, str(label)+'_full', writer)
        # advect_rk(Sol, flux, dt, elem, step%2, ud, th, mpv, node, str(label)+'_full', writer)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_advect')

        euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_ebnaexp')
        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0, writer=writer, label=str(label)+'_after_full_step')

        if writer is not None: writer.populate(str(label)+'_after_full_step','p2_half',mpv.p2_nodes_half)

        ######################################################
        # Blending : Do blending after timestep
        ######################################################
        Sol, mpv = blending.blending_after_timestep(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt, swe_to_lake, debug)

        if c1 or c2:
            print(colored("hydrostatic to nonhydrostatic conversion...", 'blue'))

            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_half_full')

            Sol = deepcopy(Sol_half_old)
            mpv = deepcopy(mpv_half_old)
            if test_hydrob == False:
                print(colored("test_hydrob == False", 'red'))
                writer.write_all(Sol,mpv,elem,node,th,str(label)+'_quarter')

                print("quarter dt = %.8f" %(dt*0.5))

                ret = time_update(Sol_half_old,flux_half_old,mpv_half_old, dt-0.5*dt, dt+0.5*dt, ud, elem, node, [0,0], th, bld=None, writer=None, debug=False)

                Sol_tu = deepcopy(ret[0])
                # mpv_tu = deepcopy(ret[2])
                Sol.rho[...] = Sol_tu.rho_half
                Sol.rhou[...] = Sol_tu.rhou_half
                Sol.rhov[...] = Sol_tu.rhov_half
                Sol.rhow[...] = Sol_tu.rhow_half
                Sol.rhoX[...] = Sol_tu.rhoX_half
                Sol.rhoY[...] = Sol_tu.rhoY_half

                # mpv.p2_nodes[...] = mpv_tu.p2_nodes_half

                writer.write_all(Sol,mpv,elem,node,th,str(label)+'_half')

                ret = time_update(Sol,flux,mpv, dt, 2.0*dt, ud, elem, node, [0,0], th, bld=None, writer=None, debug=False)

                Sol = deepcopy(ret[0])
                flux = deepcopy(ret[1])
                mpv = deepcopy(ret[2])

            if test_hydrob == True:
                # writer.write_all(Sol,mpv,elem,node,th,str(label)+'_half')

                print(colored("test_hydrob == True", 'red'))
                ret = time_update(Sol_half_old,flux_half_old,mpv_half_old, dt-0.5*dt, dt, ud, elem, node, [0,0], th, bld=None, writer=None, debug=False)

                Sol = deepcopy(ret[0])
                flux = deepcopy(ret[1])
                mpv = deepcopy(ret[2])

            # Sol_tu = deepcopy(ret[0])
            # mpv_tu = deepcopy(ret[2])
            # Sol.rho[...] = Sol_tu.rho_half
            # Sol.rhou[...] = Sol_tu.rhou_half
            # Sol.rhov[...] = Sol_tu.rhov_half
            # Sol.rhow[...] = Sol_tu.rhow_half
            # Sol.rhoX[...] = Sol_tu.rhoX_half
            # Sol.rhoY[...] = Sol_tu.rhoY_half


            # Sol.rho[...] = Sol_tu.rho_half
            # Sol.rhou[...] = Sol_tu.rhou_half
            # Sol.rhov[...] = Sol_tu.rhov_half
            # Sol.rhow[...] = Sol_tu.rhow_half
            # Sol.rhoX[...] = Sol_tu.rhoX_half
            # Sol.rhoY[...] = Sol_tu.rhoY_half

            # mpv.p2_nodes[...] = mpv_tu.p2_nodes_half

            if test_hydrob == False:
                dt *= 2.0
            if c2: ud.is_nonhydrostatic = 1

        t += dt

        if writer != None:
            writer.time = t
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')
        print("###############################################################################################")
        print("step %i done, t = %.12f, dt = %.12f, CFL = %.8f, CFL_ac = %.8f" %(step, t, dt, cfl, cfl_ac))
        print("###############################################################################################")

        step += 1
        window_step += 1


    return [Sol,flux,mpv,[window_step,step]]