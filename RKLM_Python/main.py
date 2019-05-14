import numpy as np
from management.data import data_init
from management.variable import States, Vars
from numerics_fundamentals.discretization import kgrid
from physics.gas_dynamics.thermodynamic import ThemodynamicInit
from inputs.enum_bdry import BdryType
from inputs.boundary import InitializeBdry
from physics.low_mach.mpv import MPV

from inputs.travelling_vortex_3D_48 import UserData, sol_init
from inputs.user_data import UserDataInit
from management.io import io

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

flux[0] = States(elem.sfx,ud)
if elem.ndim > 1:
    flux[1] = States(elem.sfy, ud)
if elem.ndim > 2:
    flux[2] = States(elem.sfz, ud)

th = ThemodynamicInit(ud)

# bdry = InitializeBdry(elem, ud)
mpv = MPV(elem, node, ud)

Sol0 = sol_init(Sol, mpv, elem, node, th, ud)

writer = io(ud)
writer.write_all(Sol0,mpv,elem,node,th,'000')
# Explicit_malloc
# recovery_malloc

