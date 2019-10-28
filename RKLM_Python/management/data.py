from discretization.kgrid import Grid, ElemSpaceDiscr, NodeSpaceDiscr

# dependencies of the atmospheric flow solver
from inputs.boundary import set_explicit_boundary_data
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

def data_init(ud):
    inx = ud.inx
    iny = ud.iny
    inz = ud.inz
    x0 = ud.xmin
    x1 = ud.xmax
    y0 = ud.ymin
    y1 = ud.ymax
    z0 = ud.zmin
    z1 = ud.zmax
    left = ud.bdry_type_min[0]
    right = ud.bdry_type_max[0]
    bottom = ud.bdry_type_min[1]
    top = ud.bdry_type_max[1]
    back = ud.bdry_type_min[2]
    front = ud.bdry_type_max[2]

    grid = Grid(inx,iny,inz,x0,x1,y0,y1,z0,z1,left,right,bottom,top,back,front)

    elem = ElemSpaceDiscr(grid)
    node = NodeSpaceDiscr(grid)

    return elem, node

def time_update(t,tout,ud,elem,node,step,th,Sol,flux,mpv,writer=None,debug=False):
    while ((t < tout) and (step < ud.stepmax)):
        # print("---------------------------------------")
        # print("half-time prediction of advective flux")
        # print("---------------------------------------")
        
        ud.is_compressible = is_compressible(ud,step)
        ud.is_nonhydrostatic = is_nonhydrostatic(ud,step)
        ud.nonhydrostasy = nonhydrostasy(ud,t,step)
        ud.compressibility = compressibility(ud,t,step)
        ud.acoustic_order = acoustic_order(ud,t,step)
        
        dt = dynamic_timestep(Sol,t,tout,elem,ud,th, step)

        if step < 10:
            label = '00' + str(step)
        elif step < 100:
            label = '0' + str(step)
        else:
            label = str(step)

        Sol0 = deepcopy(Sol)    
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_flux')
        
        recompute_advective_fluxes(flux, Sol)

        if debug == True: writer.populate(str(label)+'_before_advect','rhoYu',flux[0].rhoY)
        if debug == True: writer.populate(str(label)+'_before_advect','rhoYv',flux[1].rhoY)
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
            euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 1.0, label=label)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_ebnaimp')

        recompute_advective_fluxes(flux, Sol)
        mpv.p2_nodes[...] = mpv.p2_nodes0

        writer.populate(str(label)+'_after_half_step','rhoYu',flux[0].rhoY)
        writer.populate(str(label)+'_after_half_step','rhoYv',flux[1].rhoY)

        print("-----------------------------------------------")
        print("full-time step with predicted advective flux")
        print("-----------------------------------------------")

        Sol = deepcopy(Sol0)

        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_half_step')

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
            euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0)

        writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')

        t += dt
        step += 1
    return t