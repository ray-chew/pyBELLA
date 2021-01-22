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

    while ((t < tout) and (step < ud.stepmax)):
        boundary.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

        label = '%.3d' %step

        if step == 0 and writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_ic')

        dt, cfl, cfl_ac = dynamic_timestep(Sol,t,tout,elem,ud,th, step)

        if 'CFLfixed' in ud.aux:
            if step < 2 : dt = 0.2169

        ######################################################
        # Blending : Do full regime to limit regime conversion
        ######################################################
        # these make sure that we are the correct window step
        if bld is not None and window_step == 0: 
            # these make sure that blending switches are on
            if (bld.bb or bld.cb) and ud.blending_conv is not None:
                # these distinguish between SWE and Euler blending
                if ud.blending_conv == 'swe':
                    blending.do_swe_to_lake_conv(Sol, mpv, elem, node, ud, th, writer, label, debug)
                    swe_to_lake = True
                else:
                    blending.do_comp_to_psinc_conv(Sol, mpv, bld, elem, node, th, ud, label, writer)

        ######################################################
        # Blending : Do full steps or transition steps?
        ######################################################
        if bld is not None:
            c_init, c_trans = bld.criterion_init(window_step), bld.criterion_trans(window_step)
        else:
            c_init, c_trans = False, False

        ######################################################
        # Blending : If full blending steps...
        ######################################################
        # check that blending switches are on
        if c_init and bld.cb and ud.blending_conv is not None:
            # distinguish between Euler and SWE blending
            if ud.blending_conv is not 'swe':
                Sol, mpv = blending.do_psinc_to_comp_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt)
        
        ######################################################
        # Initial Blending
        ######################################################
        # Is initial blending switch on, and if yes, are we in the 0th time-step?
        if ud.initial_blending == True and step < 1 and bld is not None:
            # Distinguish between SWE and Euler blendings
            if ud.blending_conv is not 'swe':
                ud.is_compressible = 0
                ud.compressibility = 0.0
                blending.do_comp_to_psinc_conv(Sol, mpv, bld, elem, node, th, ud, label, writer)
            else:
                blending.do_swe_to_lake_conv(Sol, mpv, elem, node, ud, th, writer, label, debug)
                swe_to_lake = True
                initialise_lake_to_swe_conv = True
                ud.is_compressible = 0
                ud.compressibility = 0.0

        # Elif, is initial blending switch on and are we on the 1st time-step?
        elif ud.initial_blending == True and step == ud.no_of_pi_initial and bld is not None:
            # Distinguish between SWE and Euler blendings
            if ud.blending_conv is not 'swe':
                Sol, mpv = blending.do_psinc_to_comp_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, step, window_step, t, dt)
                ud.is_compressible = 1
                ud.compressibility = 1.0
        else:
            ud.is_compressible = is_compressible(ud,window_step)
            ud.compressibility = compressibility(ud,t,window_step)

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

        # advect(Sol, flux, 0.5*dt, elem, step%2, ud, th, mpv, node, str(label)+'_half', writer)
        advect_rk(Sol, flux, 0.5*dt, elem, step%2, ud, th, mpv, node, str(label)+'_half', writer)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_advect')

        if writer is not None: writer.populate(str(label)+'_after_full_step','p2_nodes0',mpv.p2_nodes)

        mpv.p2_nodes0[...] = mpv.p2_nodes

        euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaexp')
        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 1.0, label=str(label+'_after_ebnaimp'), writer=None)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaimp')

        recompute_advective_fluxes(flux, Sol)

        if debug == True: writer.populate(str(label)+'_after_half_step','rhoYu',flux[0].rhoY)
        if debug == True: writer.populate(str(label)+'_after_half_step','rhoYv',flux[1].rhoY)
        if debug == True and elem.ndim == 3: writer.populate(str(label)+'_after_half_step','rhoYw',flux[2].rhoY)

        mpv.p2_nodes_half = deepcopy(mpv.p2_nodes) 
        mpv.p2_nodes[...] = ud.compressibility * mpv.p2_nodes0 + (1.0-ud.compressibility) * mpv.p2_nodes
        
        Sol = deepcopy(Sol0)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_half_step')

        euler_forward_non_advective(Sol, mpv, elem, node, 0.5*dt, ud, th)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_efna')

        advect(Sol, flux, dt, elem, step%2, ud, th, mpv, node, str(label)+'_full', writer)
        # advect_rk(Sol, flux, dt, elem, step%2, ud, th, mpv, node, str(label)+'_full', writer)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_advect')

        euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_ebnaexp')
        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0, writer=None, label=str(label)+'_after_full_step')

        if writer is not None: writer.populate(str(label)+'_after_full_step','p2_half',mpv.p2_nodes_half)

        ######################################################
        # Blending : Are we in the lake regime? And is this
        #            the window step where we go back to SWE?
        ######################################################
        if bld is not None and swe_to_lake and step > 0:
            initialise_lake_to_swe_conv = bld.criterion_init(window_step+1)

        ######################################################
        # Blending : If we are in the lake regime, is blending
        #            on? If yes, do lake-to-swe conversion.
        ######################################################
        if ud.blending_conv == 'swe' and swe_to_lake and initialise_lake_to_swe_conv and bld is not None:
        # if ud.blending_conv == 'swe' and swe_to_lake and step == ud.no_of_pi_initial and bld is not None:

            # blending.do_lake_to_swe_conv(Sol, mpv, elem, node, ud, th, writer, label, debug)
            tmp_CFL = np.copy(ud.CFL)
            ud.CFL = 0.8
            Sol, mpv = blending.do_lake_to_swe_conv(Sol, flux, mpv, elem, node, ud, th, writer, label, debug, step, window_step, t, dt)
            ud.CFL = tmp_CFL
            ud.is_compressible = 1
            ud.compressibility = 1.0

        t += dt
        # print(t)

        if writer != None:
            writer.time = t
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')
        print("###############################################################################################")
        print("step %i done, t = %.12f, dt = %.12f, CFL = %.8f, CFL_ac = %.8f" %(step, t, dt, cfl, cfl_ac))
        print("###############################################################################################")

        step += 1
        window_step += 1

    return [Sol,flux,mpv,[window_step,step]]