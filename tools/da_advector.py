import sys
import os
sys.path.insert(0, os.getcwd() + '/RKLM_Python/')

from inputs import travelling_vortex_3D as tv
from inputs.user_data import UserDataInit
from management import data as data
from management import variable as var 
from management import io as io

from physics.gas_dynamics import explicit as exp
from physics.gas_dynamics import thermodynamic as thmdyn
from physics.gas_dynamics import numerical_flux as nf
from physics.low_mach import mpv as m

from data_assimilation import utils
from data_assimilation import params
from data_assimilation import letkf
from data_assimilation import etpf

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from time import time


######################################
#
# initialise arrays
#
initial_data = vars(tv.UserData())
ud = UserDataInit(**initial_data)
elem, node = data.data_init(ud)

Sol = var.Vars(elem.sc, ud)
sol_init = tv.sol_init

flux = np.empty((3), dtype=object)
flux[0] = var.States(elem.sfx, ud)
if elem.ndim > 1:
    flux[1] = var.States(elem.sfy, ud)
if elem.ndim > 2:
    flux[2] = var.States(elem.sfz, ud)

th = thmdyn.ThemodynamicInit(ud)
mpv = m.MPV(elem, node, ud)
#
######################################


######################################
#
# initialise simulation parameters
#
N = 2                   # no. of ensemble members.
step = 0
t = 0.0
dt = 0.01              # temporal step-size

label_type = 'STEP'     # output time-series in terms time or step?
#
######################################

######################################
#
# initialise data assimilation
# parameters
#
dap = params.da_params(N,da_type='etpf')
obs = dap.load_obs(dap.obs_path)
#
######################################


######################################
#
# initialise IO writer
#
ud.output_suffix = '_ensemble=%i%s' %(N, ud.output_suffix)

writer = io.io(ud)
writer.write_attrs()
if N > 1:
    writer.write_da_attrs(dap)
#
######################################


######################################
#
# initialise ensemble
#
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
ens = utils.ensemble(sol_ens)
#
######################################


######################################
#
# write initial conditions
#
for n in range(N):
    Sol = ens.members(ens)[n][0]
    mpv = ens.members(ens)[n][2]

    ens_str = 'ensemble_mem=%i_' %n if N > 1 else ''
    # set_explicit_boundary_data(Sol, elem, ud, th, mpv)
    if label_type == 'STEP':
        label = ('%s%.3d' %(ens_str,step))
    else:
        label = ('%s%.2f' %(ens_str,0.0))
    writer.write_all(Sol,mpv,elem,node,th,str(label)+'_ic')
#
######################################


######################################
#
# start simulation
#
tic = time()
tout_cnt = 0
for tout in ud.tout:

    ######################################
    #
    # start advection routine for each 
    # ensemble member
    #
    print('##############################################')
    print('Next tout = %.3f' %tout)
    print("Starting forecast...")
    results = []
    for mem in ens.members(ens):
        Sol, flux, mpv = mem[0], mem[1], mem[2]

        window_step = 0
        u,v = np.ones_like(Sol.rho), np.ones_like(Sol.rho)

        while t < tout:
            nf.recompute_advective_fluxes(flux,Sol,u=u,v=v)
            exp.advect(Sol, flux, dt, elem, window_step%2, ud, th, mpv)
        
            t += dt
            window_step += 1
        results.append([Sol,flux,mpv])
    results = np.array(results)
    #
    # end advection routine
    #
    ######################################


    ######################################
    #
    # start data assimilation step 
    #
    if N > 1 and tout in dap.da_times:

        ######################################
        #
        # algorithm:
        # localisation w/ batch observations
        #
        if dap.da_type == 'batch_obs':
            print("Starting analysis... for batch observations")
            analysis = []
            for attr in dap.obs_attributes:
                print("Assimilating %s..." %attr)

                obs_current = np.array(obs[np.argwhere(dap.da_times == tout)[0][0]][attr])

                analysis.append(letkf.da_interface(results,obs_current,dap.inflation_factor,attr,N,ud))

            analysis = np.array(analysis)

            print("Writing analysis...")
            cnt = 0
            for attr in dap.obs_attributes:
                current = analysis[cnt]
                for n in range(N):
                    setattr(results[:,dap.loc,...][n],attr,current[n])
                cnt += 1
        #
        # end localisation w/ batch
        # observations
        #
        ###################################### 
        

        ######################################
        #
        # algorithm: ETPF
        #
        elif dap.da_type == 'etpf':
            utils.ensemble_inflation(results,dap.attributes,dap.inflation_factor,N)
            results = etpf.da_interface(results,obs,dap.obs_attributes,dap.rejuvenation_factor,dap.da_times,tout,N)
        #
        # end ETPF
        #
        ###################################### 

        else:
            assert 0, "DA type not implemented: use 'rloc' or 'batch_obs'."

    ens.set_members(results)
    #
    # end data assimilation step
    #
    ######################################


    ######################################
    #
    # start writing anlysis output
    #
    print("Starting output...")
    for n in range(N):
        Sol = ens.members(ens)[n][0]
        mpv = ens.members(ens)[n][2]

        ens_str = 'ensemble_mem=%i_' %n if N > 1 else ''
        if label_type == 'STEP':
            label = ('%s%.3d' %(ens_str,step))
        else:
            label = ('%s%.2f' %(ens_str,tout))
        writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_full_step') 
    #
    # end writing analysis output
    #
    ######################################

    step += 1
    t = tout

toc = time()
print("Time taken = %.6f" %(toc-tic))

writer.close_everything()
#
# end simulation
#
######################################