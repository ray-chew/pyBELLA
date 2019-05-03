import numpy as np
from management.data import data_init
from management.variable import States
from numerics_fundamentals.discretization import kgrid
from physics.gas_dynamics.thermodynamic import ThemodynamicInit
from input.enum_bdry import BdryType
from input.boundary import InitializeBdry
from physics.low_mach.mpv import MPV

from input.travelling_vortex_3D_48 import UserData
from input.user_data import UserDataInit

step = 0
t = 0.0

initial_data = vars(UserData())
ud = UserDataInit(**initial_data)
elem, node = data_init(ud)

Sol0 = States(elem.nc, ud)
Sol = States(elem.nc, ud)
dSol = States(elem.nc, ud)
Solk = States(int(3 * ud.ncache / 2), ud)

n_aux = node.ifx
n_aux *= node.ify if node.ndim > 1 else 1
n_aux *= node.ifz if node.ndim > 2 else 1

diss = np.zeros((elem.nc))
W0 = np.zeros((n_aux))
flux = np.empty((3), dtype=object)

flux[0] = States(elem.nfx,ud)
if elem.ndim > 1:
    flux[1] = States(elem.nfy, ud)
if elem.ndim > 2:
    flux[2] = States(elem.nfz, ud) 

th = ThemodynamicInit(ud)

bdry = InitializeBdry(elem, ud)
mpv = MPV(elem, node, ud)

# Explicit_malloc
# recovery_malloc

