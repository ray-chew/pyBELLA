import numpy as np

# dependencies of the atmospheric flow solver
from management.data import data_init, time_update
from inputs.boundary import set_explicit_boundary_data
from inputs.boundary import set_ghostnodes_p2
from management.variable import States, Vars
from physics.gas_dynamics.thermodynamic import ThemodynamicInit
from physics.gas_dynamics.eos import nonhydrostasy, compressibility, is_compressible, is_nonhydrostatic
from physics.gas_dynamics.gas_dynamics import dynamic_timestep
from physics.low_mach.mpv import MPV, acoustic_order

# dependencies of the parallelisation by dask
from dask.distributed import Client, progress

# dependencies of the data assimilation module
from data_assimilation.params import da_params
from data_assimilation.utils import ensemble, sliding_window_view, ensemble_inflation, set_p2_nodes, set_rhoY_cells
from data_assimilation.letkf import da_interface, bin_func, prepare_rloc
from data_assimilation.letkf import analysis as letkf_analysis
from data_assimilation import etpf
from data_assimilation import blending
from scipy import sparse

# input file
# from inputs.baroclinic_instability_periodic import UserData, sol_init
# from inputs.travelling_vortex_2D import UserData, sol_init
from inputs.travelling_vortex_3D import UserData, sol_init
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
output_timesteps = True
label_type = 'TIME'
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
bld = blending.Blend(ud)

print("Input file is%s" %ud.output_base_name.replace('_',' '))

##########################################################
# Data Assimilation part
N = 1
dap = da_params(N,da_type='batch_obs')
# dap = da_params(N, da_type='test')
if elem.ndim == 2:
    rloc = prepare_rloc(ud, elem, node, dap, N)

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
        sol_ens[n] = [Sol0,deepcopy(flux),mpv0,[-np.inf,step]]
else:
    sol_ens = [[sol_init(Sol, mpv, elem, node, th, ud),flux,mpv,[-np.inf,step]]]
ens = ensemble(sol_ens)

##########################################################
# Load observations
# where are my observations?
if N > 1:
    obs = dap.load_obs(dap.obs_path)

# add ensemble info to filename
ud.output_suffix = '_ensemble=%i%s' %(N, ud.output_suffix)
##########################################################

if __name__ == '__main__':
    writer = io(ud)
    writer.write_attrs()
    wrtr = None
    if N > 1:
        writer.write_da_attrs(dap)
    elif output_timesteps == True:
        wrtr = writer
    for n in range(N):
        Sol = ens.members(ens)[n][0]
        mpv = ens.members(ens)[n][2]
        if label_type == 'STEP':
            label = ('ensemble_mem=%i_%.3d' %(n,step))
        else:
            label = ('ensemble_mem=%i_%.6f' %(n,0.0))
        writer.write_all(Sol,mpv,elem,node,th,str(label)+'_ic')

    client = Client(threads_per_worker=1, n_workers=1)
    tic = time()

    #### main time looping
    tout_old = -np.inf
    tout_cnt = 0
    outer_step = 0
    for tout in ud.tout:
        futures = []

        if N > 1 :
            blend = bld if tout_old in dap.da_times else None
        else:
            blend = bld

        print('##############################################')
        print('Next tout = %.3f' %tout)
        print("Starting forecast...")
        mem_cnt = 0
        for mem in ens.members(ens):
            # future = client.submit(time_update, *[mem[0],mem[1],mem[2], t, tout, ud, elem, node, mem[3], th, bld, None, False])
            if N > 1 : mem[3][0] = 0 if tout_old in dap.da_times else mem[3][0]
            if N == 1 : mem[3][0] = 0
            print("For ensemble member = %i..." %mem_cnt)
            future = time_update(mem[0],mem[1],mem[2], t, tout, ud, elem, node, mem[3], th, blend, wrtr, debug)

            futures.append(future)
            mem_cnt += 1
        results = client.gather(futures)
        results = np.array(results)

        s_res = client.scatter(results)

        #### if observations are available, do analysis...
        if N > 1 and tout in dap.da_times:
            futures = []
            if dap.da_type == 'batch_obs':
                print("Starting analysis... for batch observations")
                for attr in dap.obs_attributes:
                    print("Assimilating %s..." %attr)

                    obs_current = np.array(obs[list(dap.da_times).index(tout)][attr])

                    # future = client.submit(da_interface, *[s_res,obs_current,dap.inflation_factor,attr,N,ud,dap.loc[attr]])

                    future = da_interface(results,obs_current,dap.inflation_factor,attr,N,ud,dap.loc[attr])

                    # if attr == 'rhoY':
                    #     set_p2_nodes(future,results,N,th,node,ud)
                    # if attr == 'p2_nodes':
                    #     set_rhoY_cells(future,results,N,th,ud)

                    futures.append(future)

                # analysis = client.gather(futures)
                analysis = futures
                # analysis = np.array(analysis)

                print("Writing analysis...")
                cnt = 0
                for attr in dap.obs_attributes:
                    current = analysis[cnt]
                    for n in range(N):
                        setattr(results[:,dap.loc[attr],...][n],attr,current[n])
                    cnt += 1

            elif dap.da_type == 'test':
                print("Starting analysis... for test algorithm")
                results = rloc.analyse(results,obs,N,tout)

                # for cnt, attr in enumerate(dap.obs_attributes):
                #     current = analysis[cnt]

                #     for n in range(N):
                #         data = current[n]
                #         data = data.reshape(Nx, Ny)
                        
                #         data = np.pad(data,2,mode='constant')

                #         setattr(results[:,dap.loc[attr],...][n],attr,data)

            elif dap.da_type == 'rloc':
                print("Starting analysis... for R-localisation")
                inner = (slice(elem.igx,-elem.igx),slice(elem.igy,-elem.igy))

                Nx = elem.iicx
                Ny = elem.iicy
                attr_len = len(obs_attributes)

                tmp = np.array([getattr(results[:,loc,...][n],obs_attributes[0])[inner] for n in range(N)])
                tmp = tmp[:,np.newaxis,...]
                obs_current = np.array(obs[np.where(np.isclose(times,tout))[0][0]][obs_attributes[0]])[inner]
                
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
                    local_ens = letkf_analysis(X[n],inflation_factor,n)
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

            elif dap.da_type == 'etpf':
                ensemble_inflation(results,dap.attributes,dap.inflation_factor,N)
                results = etpf.da_interface(results,obs,dap.obs_attributes,dap.rejuvenation_factor,dap.da_times,tout,N)

            else:
                assert 0, "DA type not implemented: use 'rloc', 'batch_obs' or 'etpf'."

        for n in range(N):
            Sol = results[n][dap.loc_c]
            set_explicit_boundary_data(Sol, elem, ud, th, mpv)
            results[n][dap.loc_c] = Sol
            p2_nodes = getattr(results[n][dap.loc_n],'p2_nodes')
            set_ghostnodes_p2(p2_nodes, node, ud)
            setattr(results[n][dap.loc_n], 'p2_nodes', p2_nodes)

        ens.set_members(results, tout)

        print("Starting output...")
        for n in range(N):
            Sol = ens.members(ens)[n][0]
            mpv = ens.members(ens)[n][2]

            if label_type == 'STEP':
                # step = ens.members(ens)[0][3]
                step = outer_step
                label = ('ensemble_mem=%i_%.3d' %(n,step))
            else:
                label = ('ensemble_mem=%i_%.3f' %(n,tout))
            writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step')

        # synchronise_variables(mpv, Sol, elem, node, ud, th)
        t = tout
        tout_old = np.copy(tout)
        print('tout = %.3f' %tout)

        tout_cnt += 1
        outer_step += 1

    toc = time()
    print("Time taken = %.6f" %(toc-tic))

    writer.close_everything()