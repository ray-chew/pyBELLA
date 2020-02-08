import sys
import os
sys.path.insert(0, os.getcwd() + '/RKLM_Python/')
print(sys.path)

from inputs import travelling_vortex_3D as tv
from inputs.user_data import UserDataInit
from management import data as data
from management import variable as var 
from management import io as io

from physics.gas_dynamics import explicit as exp
from physics.gas_dynamics import thermodynamic as thmdyn
from physics.gas_dynamics import numerical_flux as nf
from physics.low_mach import mpv as m

from data_assimilation import utils as utils

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from time import time

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

N = 1
step = 0
t = 0.0
dt = 0.001

label_type = 'STEP'

ud.output_suffix = '_ensemble=%i%s' %(N, ud.output_suffix)

writer = io.io(ud)
writer.write_attrs()
if N > 1:
    writer.write_da_attrs(da_parameters)

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

# print(ens)

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

tic = time()
tout_cnt = 0
for tout in ud.tout:
    print('##############################################')
    print('Next tout = %.3f' %tout)
    print("Starting forecast...")
    for mem in ens.members(ens):
        print(mem)
        Sol, flux, mpv = mem[0], mem[1], mem[2]
        
        window_step = 0
        u,v = np.ones_like(Sol.rho), np.ones_like(Sol.rho)

        while t < tout:
            nf.recompute_advective_fluxes(flux,Sol,u=u,v=v)
            exp.advect(Sol, flux, dt, elem, window_step%2, ud, th, mpv)
        
            t += dt
            window_step += 1

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

    step += 1
    t = tout