import numpy as np

# dependencies of the atmospheric flow solver
from management.data import data_init, time_update
from management.variable import States, Vars
from physics.gas_dynamics.thermodynamic import ThemodynamicInit
from physics.gas_dynamics.eos import nonhydrostasy, compressibility, is_compressible, is_nonhydrostatic
from physics.gas_dynamics.gas_dynamics import dynamic_timestep
from physics.low_mach.mpv import MPV, acoustic_order

# dependencies of the parallelisation by dask
from dask.distributed import Client, progress

# dependencies of the data assimilation module
from data_assimilation.inputs import da_params
from data_assimilation.utils import ensemble
from data_assimilation.letkf import letkf


# input file
from inputs.travelling_vortex_3D_48 import UserData, sol_init
# from inputs.acoustic_wave_high import UserData, sol_init
# from inputs.internal_long_wave import UserData, sol_init
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

# dt_factor = 0.5 if ud.initial_impl_Euler == True else 1.0

writer = io(ud)
# writer.write_all(Sol,mpv,elem,node,th,'000_ic')

##########################################################
# Data Assimilation part
N = 20
da_parameters = da_params(N)
aprior_error_covar = da_parameters.aprior_error_covar
sampler = da_parameters.sampler_gaussian(aprior_error_covar)
# sampler = da_parameters.sampler_none()
attributes = da_parameters.attributes
ens = ensemble()
ens.initialise_members([Sol,flux,mpv],N)

# print(ens.members(ens))
ens.ensemble_spreading(sampler,attributes)

##########################################################
# Load observations
# where are my observations?
obs_path = './output_travelling_vortex_3d_48/output_travelling_vortex_3d_48_low_mach_gravity_comp_256_256.h5'
obs_file = h5py.File(obs_path)
# which attributes do I want to observe?
obs_attributes = ['rho', 'rhou', 'rhov', 'rhoY']
# when were these observations taken?
times = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

# axis 0 stores time series
obs = np.empty(len(times), dtype=object)
t_cnt = 0
for time in times:
    # how were these dataset called?
    label = '_' + str(time) + '_' + 'after_full_step'
    # axis 1 stores the attributes
    obs[t_cnt] = []
    for attribute in obs_attributes:
        obs[t_cnt].append((obs_file[str(attribute)][str(attribute) + str(label)][:]))
    t_cnt += 1
obs = np.array(obs)
obs_file.close()

# print(np.array(obs[0]).shape)
##########################################################

assert(0)

if __name__ == '__main__':
    client = Client(threads_per_worker=2, n_workers=4)    
    tic = time()
    # assert(0)

    # main time looping
    for tout in ud.tout:
        futures = []
        for mem in ens.members(ens):
            # s_ud = client.scatter(ud)
            # time_update = time_update_wrapper(t, tout, ud, elem, node, step, th, writer=writer, debug=debug)
            future = client.submit(time_update, *[mem[0],mem[1],mem[2], t, tout, ud, elem, node, step, th, None, debug])
            # future = client.submit(time_update, mem)
            # time_update(mem[0], mem[1], mem[2], t, tout, ud, elem, node, step, th, writer, debug)
            futures.append(future)
        results = client.gather(futures)
        ens.set_members(results)
        # print(results)
        # assert(0)
        
        # synchronise_variables(mpv, Sol, elem, node, ud, th)

        t = tout
        print(tout)

        for n in range(N):
            Sol = ens.members(ens)[n][0]
            mpv = ens.members(ens)[n][2]
            label = 'ensemble_mem='+str(n)+'_'+str(t)
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')

    toc = time()
    print("Time taken = %.6f" %(toc-tic))