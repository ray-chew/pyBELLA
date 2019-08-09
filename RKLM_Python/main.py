import numpy as np
from management.data import data_init
from management.variable import States, Vars
from numerics_fundamentals.discretization import kgrid
from physics.gas_dynamics.thermodynamic import ThemodynamicInit
from physics.gas_dynamics.numerical_flux import recompute_advective_fluxes
from physics.gas_dynamics.explicit import advect
from inputs.enum_bdry import BdryType
from physics.low_mach.mpv import MPV

from inputs.travelling_vortex_3D_48 import UserData, sol_init
from inputs.user_data import UserDataInit
from management.io import io

from debug import find_nearest

np.set_printoptions(precision=18)

step = 0
t = 0.0

initial_data = vars(UserData())
ud = UserDataInit(**initial_data)
elem, node = data_init(ud)

# Sol0 = States(elem.sc, ud)
Sol = Vars(elem.sc, ud)
dSol = Vars(elem.sc, ud)
# Solk = States(int(3 * ud.ncache / 2), ud)

n_aux = node.ifx
n_aux *= node.ify if node.ndim > 1 else 1
n_aux *= node.ifz if node.ndim > 2 else 1

diss = np.zeros((elem.sc))
W0 = np.zeros((n_aux))
flux = np.empty((3), dtype=object)

flux[0] = States(elem.sfx, ud)
if elem.ndim > 1:
    flux[1] = States(elem.sfy, ud)
if elem.ndim > 2:
    flux[2] = States(elem.sfz, ud)

th = ThemodynamicInit(ud)

# bdry = InitializeBdry(elem, ud)
mpv = MPV(elem, node, ud)

Sol0 = sol_init(Sol, mpv, elem, node, th, ud)

dt_factor = 0.5 if ud.initial_impl_Euler == True else 1.0

writer = io(ud)
writer.write_all(Sol0,mpv,elem,node,th,'initial')
# Explicit_malloc
# recovery_malloc

step = 0
dt = 0.0075005354646259159

# find_nearest(Sol0.rhou, 0.50000030887122227)
# print(Sol.rhou[0])
while ((t < ud.tout) and (step < ud.stepmax)):
    if step < 10:
        label = '00' + str(step)
    elif step < 100:
        label = '0' + str(step)
    else:
        label = str(step)

    print("---------------------------------------")
    print("half-time prediction of advective flux")
    print("---------------------------------------")
    
    recompute_advective_fluxes(flux, Sol0)
    writer.populate('000','rhoYu',flux[0].rhoY)
    writer.populate('000','rhoYv',flux[1].rhoY)
    
    advect(Sol0, flux, 0.5*dt, elem, step%2, ud, th, mpv)
    writer.write_all(Sol,mpv,elem,node,th,label)    
    break