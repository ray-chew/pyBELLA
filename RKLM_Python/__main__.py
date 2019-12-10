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
from data_assimilation.utils import ensemble, sliding_window_view
from data_assimilation.letkf import da_interface, bin_func
from data_assimilation.letkf import analysis as letkf_analysis
from data_assimilation import etpf
from scipy import sparse

# input file
# from inputs.baroclinic_instability_periodic import UserData, sol_init
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

# dt_factor = 0.5 if ud.initial_impl_Euler == True else 1.0

# Baroclinic IC test code snippet
# Sol = sol_init(Sol, mpv, elem, node, th, ud)
# writer = io(ud)
# writer.write_all(Sol,mpv,elem,node,th,'000_ic')

##########################################################
# Data Assimilation part
N = 30
da_parameters = da_params(N,)
da_type = da_parameters.da_type
aprior_error_covar = da_parameters.aprior_error_covar
localisation_matrix = da_parameters.localisation_matrix

sampler = da_parameters.sampler_none()
if N > 1:
    None
    # sampler = da_parameters.sampler_perturbator(5)
    # sampler = da_parameters.sampler_gaussian(aprior_error_covar)

attributes = da_parameters.attributes

print("Generating initial ensemble...")
sol_ens = np.zeros((N), dtype=object)
np.random.seed(555)
seeds = np.random.randint(10000,size=N) if N > 1 else [None]
if seeds[0] != None:
    print("Seeds used in generating initial ensemble spread = ", seeds)
    for n in range(N):
        Sol0 = deepcopy(Sol)
        mpv0 = deepcopy(mpv)
        Sol0 = sol_init(Sol0,mpv0,elem,node,th,ud, seed=seeds[n])
        sol_ens[n] = [Sol0,deepcopy(flux),mpv0,step]
else:
    sol_ens = [[sol_init(Sol, mpv, elem, node, th, ud),flux,mpv,step]]
ens = ensemble(sol_ens)

# ens = ensemble()
# ens.initialise_members([Sol,flux,mpv],N)
# ens.ensemble_spreading(ens,sampler,attributes)

# assert(0)

##########################################################
# Load observations
# where are my observations?
if N > 1:
    # obs_path = './output_travelling_vortex/output_travelling_vortex_ensemble=1_256_256_10.0.h5'
    obs_path = './output_travelling_vortex/output_travelling_vortex_ensemble=1_32_32_10.0_truth.h5'
    obs_file = h5py.File(obs_path, 'r')
    #### which attributes do I want to observe?
    obs_attributes = ['rho', 'rhou', 'rhov']
    # obs_attributes = ['rho', 'rhou', 'rhov']
    # obs_attributes = ['rhou', 'rhov']
    # obs_attributes = ['rho']
    # where in the "solutions" container are they located? 0: Sol, 1: flux, 2: mpv
    loc = 0
    #### when were these observations taken?
    # times = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    # times = [0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0]
    # times = np.linspace(0.2,10.0,99)
    # times = np.linspace(0.2,10.0,99*2-1)
    # times = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
    # print(obs_file['rho'].keys())
    # assert(0)
    steps = np.arange(0,321,8)
    times = steps / 32
    times = times[1:]

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

if elem.ndim == 2: 
    ud.output_name_comp = ("_ensemble=%i_%i_%i_%.1f" %(N,elem.icx-2*elem.igx,elem.icy-2*elem.igy,ud.tout[-1]))
if elem.ndim == 3:
    ud.output_name_comp = ("_ensemble=%i_%i_%i_%i_%.1f" %(N,elem.icx-2*elem.igx,elem.icy-2*elem.igy,elem.icz-2*elem.igz,ud.tout[-1]))
if len(ud.output_suffix) > 0:
    ud.output_name_comp += '_' + ud.output_suffix
##########################################################
# steps = np.zeros((N))
# assert(0)


if __name__ == '__main__':
    writer = io(ud)
    for n in range(N):
        Sol = ens.members(ens)[n][0]
        mpv = ens.members(ens)[n][2]
        label = ('ensemble_mem=%i_%.3d' %(n,step))
        writer.write_all(Sol,mpv,elem,node,th,str(label)+'_ic')

    client = Client(threads_per_worker=4, n_workers=2)
    tic = time()

    # assert(0)

    #### main time looping
    tout_cnt = 0
    for tout in ud.tout:
        futures = []
        print('##############################################')
        print('Next tout = %.3f' %tout)
        # obs_current = obs[tout_cnt]
        print("Starting forecast...")
        for mem in ens.members(ens):
            # s_ud = client.scatter(ud)
            # time_update = time_update_wrapper(t, tout, ud, elem, node, step, th, writer=writer, debug=debug)
            future = client.submit(time_update, *[mem[0],mem[1],mem[2], t, tout, ud, elem, node, mem[3], th, writer, debug])
            # future = client.submit(time_update, mem)
            # time_update(mem[0], mem[1], mem[2], t, tout, ud, elem, node, step, th, writer, debug)

            futures.append(future)
        results = client.gather(futures)
        results = np.array(results)
        # print(results.flatten())
        # assert(0)
        # print(results[:,loc,...])
        results_before = deepcopy(results)

        #### if observations are available, do analysis...
        # print(np.array(np.where(np.isclose(times,tout))[0]))
        # print(np.isclose(times,tout))
        # print(tout, times)
        if N > 1:
            if len(np.where(np.isclose(times,tout))[0]) > 0:
                # print(np.where(times == tout)[0][0])
                futures = []
                # anals = []
                if da_type == 'batch_obs':
                    print("Starting analysis... for batch observations")
                    for attr in obs_attributes:
                        print("Assimilating %s..." %attr)
                        obs_current = np.array(obs[np.where(np.isclose(times,tout))[0][0]][attr])
                        future = client.submit(da_interface, *[results,obs_current,attr,N,ud])
                        # analysis = da_interfaces(results,obs_current,attr,N,ud)
                        # anals.append(analysis)
                        futures.append(future)

                    analysis = client.gather(futures)
                    analysis = np.array(analysis)

                    print("Writing analysis...")
                    cnt = 0
                    for attr in obs_attributes:
                        current = analysis[cnt]
                        for n in range(N):
                            setattr(results[:,loc,...][n],attr,current[n])
                        cnt += 1

                elif da_type == 'rloc':
                    print("Starting analysis... for R-localisation")
                    inner = (slice(elem.igx,-elem.igx),slice(elem.igy,-elem.igy))

                    Nx = elem.icx - 2*elem.igx
                    Ny = elem.icy - 2*elem.igy
                    attr_len = len(obs_attributes)

                    tmp = np.array([getattr(results[:,loc,...][n],obs_attributes[0])[inner] for n in range(N)])
                    tmp = tmp[:,np.newaxis,...]
                    obs_current = np.array(obs[np.where(np.isclose(times,tout))[0][0]][obs_attributes[0]])[inner]
                    
                    # obs_current = bin_func(obs_current,(Nx,Ny))
                    obs_current = obs_current[np.newaxis,...]
                    
                    for attr in obs_attributes[1:]:
                        tmp = np.hstack((tmp,np.array([getattr(results[:,loc,...][n],attr)[inner] for n in range(N)])[:,np.newaxis,...]))
                        tmp01 = np.array(obs[np.where(np.isclose(times,tout))[0][0]][attr])[inner]
                        # tmp01 = bin_func(tmp01,(Nx,Ny))
                        # print(tmp01.shape)
                        tmp01 = tmp01[np.newaxis,...]
                        obs_current = np.vstack((obs_current,tmp01))

                    # print("obs_c (before win) = ", obs_current.shape)
                    X = np.array([sliding_window_view(mem, (1,1), (1,1)).reshape(Nx*Ny,attr_len) for mem in tmp])
                    X = np.swapaxes(X,0,1)

                    # print(tmp.shape)
                    # print(X.shape
                    
                    obs_X = 5
                    obs_Y = 5
                    #.reshape(Nx*Ny,attr_len)
                    obs_current = np.array([np.pad(obs_current_attr,2,mode='wrap') for obs_current_attr in obs_current])
                    obs_current = np.array([sliding_window_view(obs_current_attr, (obs_X,obs_Y), (1,1)) for obs_current_attr in obs_current])
                    obs_current = np.swapaxes(obs_current,0,2)
                    obs_current = np.swapaxes(obs_current,0,1)
                    obs_current = obs_current.reshape(Nx*Ny,attr_len,obs_X,obs_Y)
                    # print(obs_current.shape)
                    # assert(0)
                    
                    # print("obs_c = ", obs_current.shape)

                    tmp = np.array([getattr(results[:,loc,...][n],obs_attributes[0]) for n in range(N)])
                    tmp = tmp[:,np.newaxis,...]
                    for attr in obs_attributes[1:]:
                        tmp = np.hstack((tmp,np.array([getattr(results[:,loc,...][n],attr) for n in range(N)])[:,np.newaxis,...]))
                    Y = np.array([sliding_window_view(mem, (obs_X,obs_Y), (1,1)).reshape(Nx*Ny,attr_len,obs_X,obs_Y) for mem in tmp])
                    Y = np.swapaxes(Y,0,1)
                    # print(Y.shape)

                    obs_covar = sparse.eye(attr_len*obs_X**2,attr_len*obs_Y**2, format='csr')
                    analysis = np.empty_like(X)
                    for n in range(Nx*Ny):
                        forward_operator = lambda ensemble : Y[n]
                        local_ens = letkf_analysis(X[n],n)
                        local_ens.forward(forward_operator)
                        local_ens.localisation_matrix = localisation_matrix
                        analysis_ens = local_ens.analyse(obs_current[n],obs_covar)
                        # print(X.shape)
                        # local_ens.ensemble = np.array([np.pad(mem,2,mode='wrap') for mem in local_ens.ensemble])
                        analysis[n] = analysis_ens

                    # analysis = np.swapaxes(analysis,0,1)
                    # analysis = np.array([np.pad()])

                    analysis = np.swapaxes(analysis,0,2)
                    # print(analysis.shape)
                    cnt = 0
                    for attr in obs_attributes:
                        current = analysis[cnt]

                        for n in range(N):
                            data = current[n]
                            data = data.reshape(Nx,Ny)
                            data = np.pad(data,2,mode='wrap')
                            # print("data shape = ", data.shape)
                            
                            setattr(results[:,loc,...][n],attr,data)
                        cnt += 1

                elif da_type == 'etpf':
                    results = etpf.da_interface(results,obs,obs_attributes,times,tout,N)

                    # assert(0)
                else:
                    assert 0, "DA type not implemented: use 'rloc' or 'batch_obs'."

        # print(results_before[:,loc,...][0].rho)
        # print(results[:,loc,...][0].rho)
            # assert(0,"Assimilation failed")
        ens.set_members(results)

        # if N > 1:
        #     if np.allclose(ens.members(ens)[0].rho, results_before[:,loc,...][0].rho):
        #         print("Assimilation check: Rho quantities unchanged, i.e. no assimilation took place in rho.")

        print("Starting output...")
        for n in range(N):
            Sol = ens.members(ens)[n][0]
            mpv = ens.members(ens)[n][2]
            set_explicit_boundary_data(Sol, elem, ud, th, mpv)
            label = ('ensemble_mem=%i_%.2f' %(n,tout))
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')

        # synchronise_variables(mpv, Sol, elem, node, ud, th)
        t = tout
        print('tout = %.3f' %tout)

        tout_cnt += 1

    toc = time()
    print("Time taken = %.6f" %(toc-tic))

    writer.close_everything()