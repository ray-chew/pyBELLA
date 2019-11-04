import numpy as np

# dependencies of the atmospheric flow solver
from management.data import data_init, time_update
from inputs.boundary import set_explicit_boundary_data
from management.variable import States, Vars
from physics.gas_dynamics.thermodynamic import ThemodynamicInit
from physics.gas_dynamics.eos import nonhydrostasy, compressibility, is_compressible, is_nonhydrostatic
from physics.gas_dynamics.gas_dynamics import dynamic_timestep
from physics.low_mach.mpv import MPV, acoustic_order

# dependencies of the parallelisation by dask
from dask.distributed import Client, progress

# dependencies of the data assimilation module
from data_assimilation.params import da_params
from data_assimilation.utils import ensemble
from data_assimilation.letkf import da_interface

# input file
from inputs.travelling_vortex_3D_48 import UserData, sol_init
# from inputs.acoustic_wave_high import UserData, sol_init
# from inputs.internal_long_wave import UserData, sol_init
# from inputs.rising_bubble import UserData, sol_init
from inputs.user_data import UserDataInit
from management.io import io
import h5py

# some diagnostics and 
from copy import deepcopy
from management.debug import find_nearest
from time import time


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
# writer.write_all(Sol,mpv,elem,node,th,'000_ic')

##########################################################
# Data Assimilation part
N = 10
da_parameters = da_params(N)
aprior_error_covar = da_parameters.aprior_error_covar
# sampler = da_parameters.sampler_gaussian(aprior_error_covar)
sampler = da_parameters.sampler_none()
attributes = da_parameters.attributes
ens = ensemble()
ens.initialise_members([Sol,flux,mpv],N)
ens.ensemble_spreading(ens,sampler,attributes)

##########################################################
# Load observations
# where are my observations?
obs_path = './output_travelling_vortex/output_travelling_vortex_ensemble=1_256_256_1.0.h5'
# obs_path = './output_travelling_vortex/output_travelling_vortex_3d_48_low_mach_gravity_comp_256_256_old.h5'
obs_file = h5py.File(obs_path, 'r')
#### which attributes do I want to observe?
# obs_attributes = ['rho', 'rhou', 'rhov', 'rhoY']
# obs_attributes = ['rho', 'rhou', 'rhov']
# obs_attributes = ['rhou', 'rhov']
obs_attributes = ['rho']
# obs_attributes = ['rhoX']
# where in the "solutions" container are they located? 0: Sol, 1: flux, 2: mpv
loc = 0
#### when were these observations taken?
# times = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
times = [0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0]
# times= [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
# print(obs_file['rho'].keys())
# assert(0)

#### axis 0 stores time series
obs = np.empty(len(times), dtype=object)
t_cnt = 0
for t in times:
    #### how were these dataset called?
    label = '_ensemble_mem=0_%.2f_after_full_step' %t
    #### axis 1 stores the attributes
    obs[t_cnt] = {}
    for attribute in obs_attributes:
        dict_attr = {
            attribute: obs_file[str(attribute)][str(attribute) + str(label)][:]
        }
        obs[t_cnt].update(dict_attr)
    t_cnt += 1
obs = np.array(obs)
obs_file.close()
ud.output_name_comp = ("_ensemble=%i_%i_%i_%.1f" %(N,elem.icx-2*elem.igx,elem.icy-2*elem.igy,ud.tout[-1])) 
##########################################################

writer = io(ud)

if __name__ == '__main__':
    client = Client(threads_per_worker=4, n_workers=2)    
    tic = time()
    # assert(0)

    #### main time looping
    tout_cnt = 0
    for tout in ud.tout:
        futures = []
        print('##############################################')
        print('current tout = %.3f' %tout)
        # obs_current = obs[tout_cnt]
        print("Starting forecast...")
        for mem in ens.members(ens):
            # s_ud = client.scatter(ud)
            # time_update = time_update_wrapper(t, tout, ud, elem, node, step, th, writer=writer, debug=debug)
            future = client.submit(time_update, *[mem[0],mem[1],mem[2], t, tout, ud, elem, node, step, th, None, debug])
            # future = client.submit(time_update, mem)
            # time_update(mem[0], mem[1], mem[2], t, tout, ud, elem, node, step, th, writer, debug)
            futures.append(future)
        results = client.gather(futures)
        results = np.array(results)
        # print(results[:,loc,...])

        #### if observations are available, do analysis...
        # print(np.array(np.where(np.isclose(times,tout))[0]))
        # print(np.isclose(times,tout))
        # print(tout, times)
        if len(np.where(np.isclose(times,tout))[0]) > 0 and N > 1:
            print("Starting analysis...")
            # print(np.where(times == tout)[0][0])
            futures = []
            for attr in obs_attributes:
                obs_current = np.array(obs[np.where(np.isclose(times,tout))[0][0]][attr])
                future = client.submit(da_interface, *[results,obs_current,attr,N,ud])
                futures.append(future)

            analysis = client.gather(futures)
            analysis = np.array(analysis)
            iterable = iter(analysis.flatten())
            analysis = dict(zip(iterable, iterable))
            
            for attr in obs_attributes:
                current = analysis[attr]
                for n in range(N):
                    # results[:,loc,...][n].attr = current[n]
                    # print(attr)
                    # print(results[:,loc,...][n].rho)
                    setattr(results[:,loc,...][n],attr,current[n])
             
        ens.set_members(results)

        print("Starting output...")
        for n in range(N):
            Sol = ens.members(ens)[n][0]
            mpv = ens.members(ens)[n][2]
            # set_explicit_boundary_data(Sol, elem, ud, th, mpv)
            label = ('ensemble_mem=%i_%.2f' %(n,tout))
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')

        # synchronise_variables(mpv, Sol, elem, node, ud, th)
        t = tout
        print('tout = %.3f' %tout)

        tout_cnt += 1

    toc = time()
    print("Time taken = %.6f" %(toc-tic))

    writer.close_everything()