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
    # left = ud.bdry_type_min[0]
    # right = ud.bdry_type_max[0]
    # bottom = ud.bdry_type_min[1]
    # top = ud.bdry_type_max[1]
    # back = ud.bdry_type_min[2]
    # front = ud.bdry_type_max[2]

    # grid = Grid(inx,iny,inz,x0,x1,y0,y1,z0,z1,left,right,bottom,top,back,front)
    grid = Grid(inx,iny,inz,x0,x1,y0,y1,z0,z1)

    elem = ElemSpaceDiscr(grid)
    node = NodeSpaceDiscr(grid)

    return elem, node

# def time_update_wrapper(t,tout,ud,elem,node,step,th,writer=None,debug=False):
#     return lambda mem: time_update(mem[0],mem[1],mem[2],t,tout,ud,elem,node,step,th,writer=writer,debug=debug)

def time_update(Sol,flux,mpv,t,tout,ud,elem,node,step,th,bld=None,writer=None,debug=False):
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
        A list of `[Sol,flux,mpv]` data containers at time `tout`.
    """
    window_step, cnt_p2 = 0, 0

    while ((t < tout) and (step < ud.stepmax)):
        boundary.set_explicit_boundary_data(Sol, elem, ud, th, mpv)
        # print("---------------------------------------")
        # print("half-time prediction of advective flux")
        # print("---------------------------------------")

        label = '%.3d' %step
        dt, cfl, cfl_ac = dynamic_timestep(Sol,t,tout,elem,ud,th, step)

        if bld != None:
            c_init, c_trans = bld.criterion_init(window_step), bld.criterion_trans(window_step)
        else:
            c_init, c_trans = False, False
        
        if ud.continuous_blending == True:
            print("step = %i, window_step = %i" %(step,window_step))
            print("is_compressible = %i, is_nonhydrostatic = %i" %(ud.is_compressible, ud.is_nonhydrostatic))
            print("compressibility = %.3f, nonhydrostasy = %.3f" %(ud.compressibility,ud.nonhydrostasy))
            print("-------")

        if c_init:
            Sol_freeze = deepcopy(Sol)
            mpv_freeze = deepcopy(mpv)

            ret = time_update(Sol,flux,mpv, t, t+dt, ud, elem, node, step-1, th, bld=None, writer=None, debug=False)
            dp2n = (1.0 * ret[2].p2_nodes + 3.0 * mpv_freeze.p2_nodes) * 0.25

            # dp2n = ret[2].p2_nodes
            # ret = half_step(Sol, flux, mpv, t, t+0.5*dt, ud, elem, node, step-1, th)
            # dp2n = ret.p2_nodes

            # ret = time_update(Sol,flux,mpv, t, t+dt, ud, elem, node, step-1, th, bld=None, writer=None, debug=False)
            # dp2n = (1.0 * ret[2].p2_nodes + 1.0 * mpv_freeze.p2_nodes) * 0.5


            Sol = Sol_freeze
            mpv = mpv_freeze

            print("Blending... step = %i" %window_step)
            writer.populate(str(label)+'_before_blending', 'dp2n', dp2n)
            bld.convert_p2n(dp2n)
            bld.update_Sol(Sol,elem,node,th,ud, mpv,label=label,writer=writer)
            bld.update_p2n(Sol,mpv,node,th,ud)
        
        # if step > 0 and bld != None:
        #     Sol_freeze = deepcopy(Sol)
        #     mpv_freeze = deepcopy(mpv)
        #     step_tmp = np.copy(step)
        #     ret = time_update(Sol,flux,mpv, t, t+0.5*dt, ud, elem, node, step-1, th, bld=None, writer=None, debug=False)
        #     dp2n = (ret[2].p2_nodes + mpv_freeze.p2_nodes) * 0.5
        #     # dp2n = ret[2].p2_nodes

        #     Sol = Sol_freeze
        #     mpv = mpv_freeze
        #     step = step_tmp

        ud.is_compressible = is_compressible(ud,step)
        ud.is_nonhydrostatic = is_nonhydrostatic(ud,step)
        ud.nonhydrostasy = nonhydrostasy(ud,t,step)
        ud.compressibility = compressibility(ud,t,step)
        ud.acoustic_order = acoustic_order(ud,t, window_step)

        if step == 0 and writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_ic')

        Sol0 = deepcopy(Sol)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_flux')
        
        recompute_advective_fluxes(flux, Sol)

        if debug == True: writer.populate(str(label)+'_before_advect','rhoYu',flux[0].rhoY.T)
        if debug == True: writer.populate(str(label)+'_before_advect','rhoYv',flux[1].rhoY.T)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_advect')

        advect(Sol, flux, 0.5*dt, elem, step%2, ud, th, mpv)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_advect')

        mpv.p2_nodes0[...] = mpv.p2_nodes

        if ud.is_ArakawaKonor:
            ud.is_nonhydrostatic = 0
            ud.nonhydrostasy = 0.0
            ud.is_compressible = 1
            ud.compressibility = 1.0

            Sol_tmp = deepcopy(Sol)
            euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
            if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaexp')
            euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 1.0, label=label)

            ud.is_nonhydrostatic = 1
            ud.nonhydrostasy = 1.0
            ud.is_compressible = 0
            ud.compressibility = 0.0

            Sol = Sol_tmp
            euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
            if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaexp')
            euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 1.0, label=label)

        else:
            euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
            if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaexp')
            euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 1.0, label=str(label+'_after_ebnaimp'), writer=writer)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaimp')

        recompute_advective_fluxes(flux, Sol)

        if debug == True: writer.populate(str(label)+'_after_half_step','rhoYu',flux[0].rhoY.T)
        if debug == True: writer.populate(str(label)+'_after_half_step','rhoYv',flux[1].rhoY.T)

        p2_half = np.copy(mpv.p2_nodes)

        mpv.p2_nodes[...] = ud.compressibility * mpv.p2_nodes0 + (1.0-ud.compressibility) * mpv.p2_nodes

        Sol = deepcopy(Sol0)

        # print("step = %i, compressibility = %.4f, is_compressible = %i" %(window_step, ud.compressibility, ud.is_compressible))

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_half_step')

        # print("-----------------------------------------------")
        # print("full-time step with predicted advective flux")
        # print("-----------------------------------------------")

        euler_forward_non_advective(Sol, mpv, elem, node, 0.5*dt, ud, th)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_efna')

        advect(Sol, flux, dt, elem, step%2, ud, th, mpv)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_advect')

        if ud.is_ArakawaKonor:
            ud.is_nonhydrostatic = 0
            ud.nonhydrostasy = 0.0
            ud.is_compressible = 1
            ud.compressibility = 1.0

            Sol_tmp = deepcopy(Sol)
            euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
            euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0)

            ud.is_nonhydrostatic = 1
            ud.nonhydrostasy = 1.0
            ud.is_compressible = 0
            ud.compressibility = 0.0

            Sol = Sol_tmp
            euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
            if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_ebnaexp')
            euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0)

        else:
            euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
            if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_ebnaexp')
            euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0, writer=writer, label=str(label)+'_after_full_step')

        # tmp = np.copy(mpv.p2_nodes)
        # if step > 0 and bld != None:
        # mpv.p2_nodes[...] = dp2n
        # mpv.p2_nodes[...] = p2_half
        if writer != None:
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')
            print("###############################################################################################")
            print("step %i done, t = %.12f, dt = %.12f, CFL = %.8f, CFL_ac = %.8f" %(step, t, dt, cfl, cfl_ac))
            print("###############################################################################################")
        t += dt
        step += 1
        window_step += 1
        # mpv.p2_nodes[...] = tmp

        # print(cnt_p2)
        # print(window_step)
        # print(t, step)
        # if step == 10:
        #     assert(0)

    return [Sol,flux,mpv,step]

def half_step(Sol,flux,mpv,t,tout,ud,elem,node,step,th):
    dt, _, _ = dynamic_timestep(Sol,t,tout,elem,ud,th, step)
    recompute_advective_fluxes(flux, Sol)
    advect(Sol, flux, 0.5*dt, elem, step%2, ud, th, mpv)

    euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
    euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 1.0, label=None, writer=None)

    return mpv
