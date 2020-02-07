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

def time_update(Sol,flux,mpv,t,tout,ud,elem,node,step,th,writer=None,debug=False):
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
    window_step = 0

    while ((t < tout) and (step < ud.stepmax)):
        boundary.set_explicit_boundary_data(Sol, elem, ud, th, mpv)
        # print("---------------------------------------")
        # print("half-time prediction of advective flux")
        # print("---------------------------------------")
        # assert(0)
        ud.is_compressible = is_compressible(ud,window_step)
        ud.is_nonhydrostatic = is_nonhydrostatic(ud,window_step)
        ud.nonhydrostasy = nonhydrostasy(ud,t,window_step)
        ud.compressibility = compressibility(ud,t,window_step)
        ud.acoustic_order = acoustic_order(ud,t,window_step)

        if ud.continuous_blending == True:
            print("is_compressible = %i, is_nonhydrostatic = %i" %(ud.is_compressible, ud.is_nonhydrostatic))
            print("nonhydrostasy = %.3f, compressibility = %.3f" %(ud.nonhydrostasy, ud.compressibility))

        dt, cfl, cfl_ac = dynamic_timestep(Sol,t,tout,elem,ud,th, step)

        label = '%.3d' %step
        
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

        # print(ud.Msq)
        # Y = Sol.rhoY / Sol.rho
        # print(Y)
        # print(mpv.HydroState_n.rhoY0)

        # offset = 0 if step % 2 == 0 and step == (ud.no_of_pi_initial + ud.no_of_pi_transition) and step != 0 else 0
        offset = 0

        c1 = step == (ud.no_of_pi_initial) and ud.continuous_blending == True and ud.no_of_pi_initial > 0

        c2 = step <= (ud.no_of_pi_transition) and ud.continuous_blending == True and ud.no_of_pi_transition > 0

        if c1 or c2:
        # if False:
            # assert(0)
            print("WINDOW STEP === %i, STEP == %i" %(window_step,step))
            print("OFFSET == %i" %offset)
            # Sol = deepcopy(Sol0)
            # mpv.p2_nodes[...] = mpv.p2_nodes0
            # mpv.p2_nodes[...] = mpv.p2_nodes
            ndim = elem.ndim
            # p2 = 0.5*(mpv.p2_nodes0 - mpv.p2_nodes0.mean()) + 0.5*(mpv.p2_nodes - mpv.p2_nodes.mean())
            if c1:
                
                p2 = 0.5*(mpv.p2_nodes0 + mpv.p2_nodes)
            # p2 = mpv.p2_nodes
                dp2n_psinc = p2 #- p2.mean()
            else:
                dp2n_psinc = mpv.p2_nodes0 #- mpv.p2_nodes0
            # dp2n_psinc = np.copy(mpv.p2_nodes)
            # dp2n_psinc = mpv.p2_nodes0 - mpv.p2_nodes0.mean()
            # dp2n_psinc = mpv.p2_nodes0 - mpv.p2_nodes0.mean()

            # mpv.p2_nodes[...] = dp2n_psinc

            kernel = np.ones([2]*ndim)
            dp2c_psinc = signal.fftconvolve(dp2n_psinc,kernel, mode='valid') / kernel.sum()
            # writer.populate(str(label)+'_after_full_step','dp2n_psinc',dp2n_psinc)
            # writer.populate(str(label)+'_after_full_step','dp2c_psinc',dp2c_psinc)

            # rhoYn = np.zeros_like(dp2n_psinc)
            # i1 = [slice(None,)]*ndim
            # igs = node.igs
            # for dim in range(ndim):
            #     i1[dim] = slice(igs[dim]-1,-igs[dim]+1)
            # i1 = tuple(i1)
            # rhoYn[i1] = signal.fftconvolve(Sol.rhoY,kernel,mode='valid') / kernel.sum()
            # boundary.set_ghostnodes_p2(rhoYn,node,ud)
            # writer.populate(str(label)+'_after_full_step','rhoYc',Sol.rhoY)
            # writer.populate(str(label)+'_after_full_step','rhoYn',rhoYn)
            # print((rhoYn**th.gm1).mean())
            # print(ud.Msq)
            # Sol = deepcopy(Sol0)
            

            # euler_forward_non_advective(Sol, mpv, elem, node, 0.5*dt, ud, th)

            # advect(Sol, flux, dt, elem, step%2, ud, th, mpv)

            # euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
            # euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0, writer=None, label=str(label)+'_after_full_step')
            if c1:
                fac = ud.Msq
                # mpv.p2_nodes[...] =  1./fac * (rhoYn**th.gm1 - 1.0 + fac * dp2n_psinc)
                mpv.p2_nodes[...] = dp2n_psinc
                Y = Sol.rhoY / Sol.rho
                # # print(Y)

                # Yn = np.zeros_like(dp2n_psinc)
                # Yn[i1] = signal.fftconvolve(Y,kernel,mode='valid') / kernel.sum()
                # boundary.set_ghostnodes_p2(Yn,node,ud)
                
                Sol = deepcopy(Sol0)

                # writer.populate(str(label)+'_after_full_step','Yn',Yn)
                rho_psinc = np.copy(Sol.rho)

                Sol.rhoY[...] = (Sol.rhoY**th.gm1 + fac * dp2c_psinc)**(th.gm1inv)
                Sol.rho[...] = Sol.rhoY / Y

                rho_fac = Sol.rho / rho_psinc
                Sol.rhou[...] *= rho_fac
                Sol.rhov[...] *= rho_fac
                Sol.rhow[...] *= rho_fac
                Sol.rhoX[...] *= rho_fac
                boundary.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

            # rhoYncomp = np.zeros_like(dp2n_psinc)
            # rhoYncomp[i1] = signal.fftconvolve(Sol.rhoY,kernel,mode='valid') / kernel.sum()
            # boundary.set_ghostnodes_p2(rhoYncomp,node,ud)
            # # mpv.p2_nodes[...] = (rhoYncomp**th.gm1 - 1.0)/fac + Yn
            # writer.populate(str(label)+'_after_full_step','rhoYcomp',rhoYncomp)
        # mpv.p2_nodes0 = mpv.p2_nodes# - mpv.p2_nodes.mean()

        # else:
        # if step == (ud.no_of_pi_initial + ud.no_of_pi_transition) and ud.continuous_blending == True:
        # if step >= 0:
        #     mpv.p2_nodes[...] = mpv.p2_nodes0
        # elif step == 0:
        #     mpv.p2_nodes[...] = 0.5 * (mpv.p2_nodes0 + mpv.p2_nodes)
            # mpv.p2_nodes -= mpv.p2_nodes.mean()
        # elif step == 0:
        #     mpv.p2_nodes[...] = mpv.p2_nodes

        # elif step <= ud.no_of_pi_transition and ud.no_of_pi_transition > 0:
        #     mpv.p2_nodes[...] = ud.compressibility * mpv.p2_nodes0 + (1.0-ud.compressibility) * mpv.p2_nodes
            
        #     dp2n_psinc = np.copy(mpv.p2_nodes)
        #     kernel = np.ones([2]*dp2n_psinc.ndim)
        #     dp2c_psinc = signal.fftconvolve(dp2n_psinc,kernel, mode='valid') / kernel.sum()

        #     Y = Sol.rhoY / Sol.rho
        #     Sol = deepcopy(Sol0)

        #     rho_psinc = np.copy(Sol.rho)

        #     Sol.rhoY[...] = (Sol.rhoY**th.gm1 + ud.Msq * dp2c_psinc)**(th.gm1inv)
        #     Sol.rho[...] = Sol.rhoY / Y

        #     rho_fac = Sol.rho / rho_psinc
        #     Sol.rhou[...] *= rho_fac
        #     Sol.rhov[...] *= rho_fac
        #     Sol.rhow[...] *= rho_fac
        #     Sol.rhoX[...] *= rho_fac
        #     boundary.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

        else:
            # if ud.is_compressible == 0 and step >= ud.no_of_pi_transition:
            #     mpv.p2_nodes[...] = mpv.p2_nodes
            # else:
                # mpv.p2_nodes[...] = mpv.p2_nodes0
                # mpv.p2_nodes[...] = 0.5 * (mpv.p2_nodes0 + mpv.p2_nodes)
            mpv.p2_nodes[...] = ud.compressibility * mpv.p2_nodes0 + (1.0-ud.compressibility) * mpv.p2_nodes
            Sol = deepcopy(Sol0)
            # mpv.p2_nodes -= mpv.p2_nodes.mean()
            # None
            # None 
            # mpv.p2_nodes[...] = mpv.p2_nodes
        # if step == 0:

        #     mpv.p2_nodes[...] = 0.5 * mpv.p2_nodes + 0.5 * mpv.p2_nodes0
        # else:
            # mpv.p2_nodes[...] = mpv.p2_nodes0
        # Sol = Sol_tmp

        print("step = %i, compressibility = %.4f, is_compressible = %i" %(window_step, ud.compressibility, ud.is_compressible))

        if step == 41:
            assert(0)

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

        if writer != None:
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')
        print("###############################################################################################")
        print("step %i done, t = %.12f, dt = %.12f, CFL = %.8f, CFL_ac = %.8f" %(step, t, dt, cfl, cfl_ac))
        print("###############################################################################################")
        t += dt
        step += 1
        window_step += 1
        # print(window_step)
        # print(t, step)
        # if step == 10:
        #     assert(0)

    return [Sol,flux,mpv,step]