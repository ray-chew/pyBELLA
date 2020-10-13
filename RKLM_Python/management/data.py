from discretization.kgrid import Grid, ElemSpaceDiscr, NodeSpaceDiscr

# dependencies of the atmospheric flow solver
import inputs.boundary as boundary
from management.variable import States, Vars
from discretization import kgrid
from physics.gas_dynamics.thermodynamic import ThemodynamicInit
from physics.gas_dynamics.numerical_flux import recompute_advective_fluxes
from physics.gas_dynamics.explicit import advect
from physics.gas_dynamics.eos import nonhydrostasy, compressibility, synchronise_variables, is_compressible, is_nonhydrostatic
from physics.gas_dynamics.gas_dynamics import dynamic_timestep
from physics.low_mach.second_projection import euler_backward_non_advective_impl_part, euler_backward_non_advective_expl_part, euler_forward_non_advective
from inputs.enum_bdry import BdryType
from physics.low_mach.mpv import MPV, acoustic_order

import numpy as np
from copy import deepcopy
from scipy import signal
import h5py

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
    th : :class:`physics.gas_dynamics.thermodynamic.ThemodynamicInit`
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

        if bld is not None and window_step == 0:
            if (bld.bb or bld.cb) and ud.blending_conv is not None:
                if ud.blending_conv == 'swe':
                    print("swe to lake conversion...")

                    setattr(ud,'min_val',Sol.rho.min()) # save the min value to Sol container
                    # eps = Sol.rho.max() - Sol.rho.min()
                    # setattr(ud,'min_val',1.0)
                    Sol.rhou[...] = Sol.rhou / Sol.rho * ud.min_val
                    # Sol.rhov[...] = Sol.rhov / Sol.rho * ud.min_val
                    Sol.rhow[...] = Sol.rhow / Sol.rho * ud.min_val
                    Sol.rhoY[...] = Sol.rhoY / Sol.rho * ud.min_val
                    Sol.rho[...] = ud.min_val

                    setattr(ud,'p2_min_val',Sol.rhoY.min())
                    # setattr(ud,'p2_min_val',np.sqrt(9.81/2.0))
                    eps = mpv.p2_nodes.max() - ud.p2_min_val
                    mpv.p2_nodes[:,1:,:] = (mpv.p2_nodes[:,1:,:] - ud.p2_min_val) #/ eps
                    # mpv.p2_nodes[:,1:,:] -= mpv.p2_nodes[:,1:,:].mean(axis=(0,2),keepdims=True)

                    if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_swe_to_lake')
                    swe_to_lake = True

                else:
                    dp2n = mpv.p2_nodes
                    bld.convert_p2n(dp2n)
                    bld.update_Sol(Sol,elem,node,th,ud,mpv,'bef',label=label,writer=writer)
                    bld.update_p2n(Sol,mpv,node,th,ud)

        if bld is not None:
            c_init, c_trans = bld.criterion_init(window_step), bld.criterion_trans(window_step)
        else:
            c_init, c_trans = False, False

        dt, cfl, cfl_ac = dynamic_timestep(Sol,t,tout,elem,ud,th, step)

        if c_init and bld.cb and ud.blending_conv is not None:
            if ud.blending_conv == 'swe':
                None
                # print("lake to swe conversion...")

                # H10 = mpv.p2_nodes - mpv.p2_nodes.min()
                # H1 = H10[:,1,:] # project horizontal slice to 2D

                # # define 2D kernel
                # kernel = np.ones((2,2))
                # kernel /= kernel.sum()

                # # do node-to-cell averaging
                # H1 = signal.convolve(H1, kernel, mode='valid')

                # # project H1 back to horizontal slice with ghost cells
                # H1 = np.expand_dims(H1, 1)
                # H1 = np.repeat(H1, elem.icy, axis=1)

                # Sol.rho += H1 #+ ud.min_val
                # Sol.rhou = Sol.rhou / ud.min_val * Sol.rho
                # # Sol.rhov = Sol.rhov / ud.min_val * Sol.rho
                # Sol.rhow = Sol.rhow / ud.min_val * Sol.rho
                # Sol.rhoY = Sol.rhoY / ud.min_val * (ud.min_val+H1)#* Sol.rho
                # mpv.p2_nodes = H10 + ud.p2_min_val

                # if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_lake_to_swe')

            else:
                print("Blending... step = %i" %window_step)
                Sol_freeze = deepcopy(Sol)
                mpv_freeze = deepcopy(mpv)

                ret = time_update(Sol,flux,mpv, t, t+dt, ud, elem, node, [0,step-1], th, bld=None, writer=None, debug=False)

                fac_old = ud.blending_weight
                fac_new = 1.0 - fac_old
                dp2n = (fac_new * ret[2].p2_nodes + fac_old * mpv_freeze.p2_nodes)
                # dp2n = mpv.p2_nodes

                if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_start', mpv_freeze.p2_nodes)
                if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_end', ret[2].p2_nodes)
                Sol = Sol_freeze
                mpv = mpv_freeze

                if writer != None: writer.populate(str(label)+'_after_full_step', 'dp2n', dp2n)
                bld.convert_p2n(dp2n)
                bld.update_Sol(Sol,elem,node,th,ud,mpv,'aft',label=label,writer=writer)
                bld.update_p2n(Sol,mpv,node,th,ud)

        if ud.initial_blending == True and step < 1 and bld is not None:
            print("passed")
            ud.is_compressible = 0
            ud.compressibility = 0.0
            dp2n = mpv.p2_nodes
            bld.convert_p2n(dp2n)
            bld.update_Sol(Sol,elem,node,th,ud,mpv,'bef',label=label,writer=writer)
            bld.update_p2n(Sol,mpv,node,th,ud)

        elif ud.initial_blending == True and step == 1 and bld is not None:
            # dp2n = mpv.p2_nodes

            Sol_freeze = deepcopy(Sol)
            mpv_freeze = deepcopy(mpv)

            ret = time_update(Sol,flux,mpv, t, t+dt, ud, elem, node, [0,step-1], th, bld=None, writer=None, debug=False)

            fac_old = ud.blending_weight
            fac_new = 1.0 - fac_old
            dp2n = (fac_new * ret[2].p2_nodes + fac_old * mpv_freeze.p2_nodes)
            # dp2n = mpv.p2_nodes

            if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_start', mpv_freeze.p2_nodes)
            if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_end', ret[2].p2_nodes)
            Sol = Sol_freeze
            mpv = mpv_freeze

            bld.convert_p2n(dp2n)
            bld.update_Sol(Sol,elem,node,th,ud,mpv,'aft',label=label,writer=writer)
            bld.update_p2n(Sol,mpv,node,th,ud)
            ud.is_compressible = 1
            ud.compressibility = 1.0
        else:
            # if ud.initial_blending == False and window_step == 0:
                # window_step = -np.inf
            ud.is_compressible = is_compressible(ud,window_step)
            ud.compressibility = compressibility(ud,t,window_step)

        ud.is_nonhydrostatic = is_nonhydrostatic(ud,window_step)
        ud.nonhydrostasy = nonhydrostasy(ud,t,window_step)
        ud.acoustic_order = acoustic_order(ud,t, window_step)

        if ud.continuous_blending == True or ud.initial_blending == True:
            if window_step >= 0:
                print("step = %i, window_step = %i" %(step,window_step))
            else:
                print("step = %i, window_step = %f" %(step,window_step))
        print("is_compressible = %i, is_nonhydrostatic = %i" %(ud.is_compressible, ud.is_nonhydrostatic))
        print("compressibility = %.3f, nonhydrostasy = %.3f" %(ud.compressibility,ud.nonhydrostasy))
        print("-------")

        if step == 0 and writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_ic')

        Sol0 = deepcopy(Sol)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_flux')
        
        recompute_advective_fluxes(flux, Sol)

        if debug == True: writer.populate(str(label)+'_before_advect','rhoYu',flux[0].rhoY)
        if debug == True: writer.populate(str(label)+'_before_advect','rhoYv',flux[1].rhoY)
        if debug == True and elem.ndim == 3: writer.populate(str(label)+'_before_advect','rhoYw',flux[2].rhoY)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_advect')

        advect(Sol, flux, 0.5*dt, elem, step%2, ud, th, mpv, node, str(label)+'_half', writer)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_advect')

        mpv.p2_nodes0[...] = mpv.p2_nodes

        euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaexp')
        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 1.0, label=str(label+'_after_ebnaimp'), writer=writer)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaimp')

        recompute_advective_fluxes(flux, Sol)

        if debug == True: writer.populate(str(label)+'_after_half_step','rhoYu',flux[0].rhoY)
        if debug == True: writer.populate(str(label)+'_after_half_step','rhoYv',flux[1].rhoY)
        if debug == True and elem.ndim == 3: writer.populate(str(label)+'_after_half_step','rhoYw',flux[2].rhoY)

        mpv.p2_nodes[...] = ud.compressibility * mpv.p2_nodes0 + (1.0-ud.compressibility) * mpv.p2_nodes
        # mpv.p2_nodes[...] = mpv.p2_nodes0

        Sol = deepcopy(Sol0)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_half_step')

        euler_forward_non_advective(Sol, mpv, elem, node, 0.5*dt, ud, th)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_efna')

        advect(Sol, flux, dt, elem, step%2, ud, th, mpv, node, str(label)+'_full', writer)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_advect')

        euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_ebnaexp')
        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0, writer=writer, label=str(label)+'_after_full_step')

        t += dt

        if writer != None:
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')
            print("###############################################################################################")
            print("step %i done, t = %.12f, dt = %.12f, CFL = %.8f, CFL_ac = %.8f" %(step, t, dt, cfl, cfl_ac))
            print("###############################################################################################")

        if bld is not None:
            lake_to_swe_init = bld.criterion_init(window_step+1)

        if ud.blending_conv == 'swe' and swe_to_lake and lake_to_swe_init:
            print("lake to swe conversion...")

            H10 = mpv.p2_nodes - mpv.p2_nodes.min()
            # H10 -= H10.mean(axis=(0,2),keepdims=True)
            H1 = H10[:,1,:] # project horizontal slice to 2D
            # H1 = H1.mean()

            # define 2D kernel
            kernel = np.ones((2,2))
            kernel /= kernel.sum()

            # do node-to-cell averaging
            H1 = signal.convolve(H1, kernel, mode='valid')

            # project H1 back to horizontal slice with ghost cells
            H1 = np.expand_dims(H1, 1)
            H1 = np.repeat(H1, elem.icy, axis=1)

            Sol.rho += H1 #+ ud.min_val
            Sol.rhou = Sol.rhou / ud.min_val * Sol.rho
            # Sol.rhov = Sol.rhov / ud.min_val * Sol.rho
            Sol.rhow = Sol.rhow / ud.min_val * Sol.rho
            Sol.rhoY = Sol.rhoY / ud.min_val * (ud.min_val+H1)#* Sol.rho
            mpv.p2_nodes = H10 + ud.p2_min_val
            # mpv.p2_nodes -= mpv.p2_nodes.mean()

            if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_lake_to_swe')

        step += 1
        # print("step = ", step)
        window_step += 1

    return [Sol,flux,mpv,[window_step,step]]