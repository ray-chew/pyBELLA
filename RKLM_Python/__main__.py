import numpy as np

from inputs.boundary import set_explicit_boundary_data
from management.data import data_init
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

# from inputs.travelling_vortex_3D_48 import UserData, sol_init
# from inputs.acoustic_wave_high import UserData, sol_init
from inputs.internal_long_wave import UserData, sol_init
# from inputs.rising_bubble import UserData, sol_init
from inputs.user_data import UserDataInit
from management.io import io
from copy import deepcopy

from management.debug import find_nearest
from time import time
import h5py

debug = False

np.set_printoptions(precision=18)

step = 0
t = 0.0

initial_data = vars(UserData())
ud = UserDataInit(**initial_data)

elem, node = data_init(ud)

Sol = Vars(elem.sc, ud)

flux = np.empty((3), dtype=object)
flux[0] = States(elem.sfx, ud)
if elem.ndim > 1:
    flux[1] = States(elem.sfy, ud)
if elem.ndim > 2:
    flux[2] = States(elem.sfz, ud)

th = ThemodynamicInit(ud)

mpv = MPV(elem, node, ud)

Sol0 = sol_init(Sol, mpv, elem, node, th, ud)
Sol = deepcopy(Sol0)

dt_factor = 0.5 if ud.initial_impl_Euler == True else 1.0

writer = io(ud)
writer.write_all(Sol,mpv,elem,node,th,'000_ic')

step = 0
tout = ud.tout

tic = time()

while ((t < tout) and (step < ud.stepmax)):
    
    dt = dynamic_timestep(Sol,t,tout,elem,ud,th, step)
    Sol0 = deepcopy(Sol)

    if step < 10:
        label = '00' + str(step)
    elif step < 100:
        label = '0' + str(step)
    else:
        label = str(step)
    
    print("---------------------------------------")
    print("half-time prediction of advective flux")
    print("---------------------------------------")
    
    ud.is_compressible = is_compressible(ud,step)
    ud.is_nonhydrostatic = is_nonhydrostatic(ud,step)
    ud.nonhydrostasy = nonhydrostasy(ud,t,step)
    ud.compressibility = compressibility(ud,t,step)
    ud.acoustic_order = acoustic_order(ud,t,step)

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

        # mpv.p2_nodesh[...] = mpv.p2_nodes
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

    print("-----------------------------------------------")
    print("full-time step with predicted advective flux")
    print("-----------------------------------------------")

    writer.populate(str(label)+'_after_half_step','rhoYu',flux[0].rhoY)
    writer.populate(str(label)+'_after_half_step','rhoYv',flux[1].rhoY)

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

        # mpv.p2_nodesh[...] = mpv.p2_nodes
        Sol = Sol_tmp
        euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_ebnaexp')
        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0)

    else:
        euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_ebnaexp')
        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0)

    writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')

    # synchronise_variables(mpv, Sol, elem, node, ud, th)
    t += dt
    step += 1
    dt_factor = 1.0

    print("############################################################################################")
    print("step %i done, t = %.12f, dt = %.12f" %(step, t, dt))
    print("############################################################################################")

toc = time()
print("Time taken = %.6f" %(toc-tic))