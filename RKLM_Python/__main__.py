import numpy as np

# dependencies of the atmospheric flow solver
from management.data import data_init, time_update, time_update_wrapper
from management.variable import States, Vars
from physics.gas_dynamics.thermodynamic import ThemodynamicInit
from physics.gas_dynamics.eos import nonhydrostasy, compressibility, is_compressible, is_nonhydrostatic
from physics.gas_dynamics.gas_dynamics import dynamic_timestep
from physics.low_mach.mpv import MPV, acoustic_order

# dependencies of the parallelisation by dask
from dask.distributed import Client, progress
# client = Client(threads_per_worker=2, n_workers=4)

# dependencies of the data assimilation module
from data_assimilation.inputs import da_params
from data_assimilation.utils import ensemble
from data_assimilation.letkf import letkf


# input file
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
# writer.write_all(Sol,mpv,elem,node,th,'000_ic')

tout = ud.tout

if __name__ == '__main__':
    ##########################################################
    # Data Assimilation part
    N = 2
    da_parameters = da_params(N)
    aprior_error_covar = da_parameters.aprior_error_covar
    # sampler = da_parameters.sampler(aprior_error_covar)
    sampler = da_parameters.sampler_none()
    attributes = da_parameters.attributes
    ens = ensemble()
    ens.initialise_members([Sol,flux,mpv],N)

    # print(ens.members(ens))
    # ens.ensemble_spreading(sampler,attributes)

    client = Client(threads_per_worker=2, n_workers=2)
    ##########################################################

    tic = time()
    # assert(0)

    # main time looping
    for tout in [ud.tout]:
        futures = []
        for mem in ens.members(ens):
            # s_ud = client.scatter(ud)
            # time_update = time_update_wrapper(t, tout, ud, elem, node, step, th, writer=writer, debug=debug)
            future = client.submit(time_update, *[mem[0], mem[1], mem[2], t, tout, ud, elem, node, step, th, None, debug])
            # future = client.submit(time_update, mem)
            # time_update(t, tout, ud, elem, node, step, th, Sol, flux, mpv, writer, debug)
            futures.append(future)
        client.gather(futures)
        # synchronise_variables(mpv, Sol, elem, node, ud, th)

        # print("############################################################################################")
        # print("step %i done, t = %.12f, dt = %.12f" %(step, t, dt))
        # print("############################################################################################")
        t = tout


    toc = time()
    print("Time taken = %.6f" %(toc-tic))