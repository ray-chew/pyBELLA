import numpy as np
from management.data import data_init
from management.variable import States
from numerics_fundamentals.discretization import kgrid
from physics.gas_dynamics.thermodynamic import ThemodynamicInit
from inputs.enum_bdry import BdryType
from inputs.boundary import InitializeBdry
from physics.low_mach.mpv import MPV

from inputs.acoustic_wave_high import UserData, sol_init
from inputs.user_data import UserDataInit

step = 0
t = 0.0

initial_data = vars(UserData())
ud = UserDataInit(**initial_data)
elem, node = data_init(ud)

# Sol0 = States(elem.sc, ud)
Sol = States(elem.sc, ud)
dSol = States(elem.sc, ud)
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

bdry = InitializeBdry(elem, ud)
mpv = MPV(elem, node, ud)

Sol = sol_init(Sol, mpv,bdry,elem, node, th, ud)
# Explicit_malloc
# recovery_malloc

