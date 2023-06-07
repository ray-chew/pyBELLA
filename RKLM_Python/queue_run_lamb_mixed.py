from run import run_params as rp
from inputs.enum_bdry import BdryType

import json
import numpy as np

rp = rp()
path_to_obs = './'

rp.N = 1

rp.tc = 'lw_p'
# rp.tc = 'mark'

t_ref = 100.0
omega = 7.292 * 1e-5
# resol_x = [151,301,601,1201]
# resol_y = [15,30,60,120]
resol_x = [301]
resol_y = [120]
resol_t = [1,2,4,8,10,12,14,16]
# resol_t = [200,400,600,800]
resol_t = [10.0]
# omegas = [0.0, 2.0 * omega * t_ref]
# omegas = [2.0 * omega * t_ref]
omegas = [0.0]


ud = {}
dap = {
    'None' : None,
}

for x,y in zip(resol_x,resol_y):
    ud['inx'] = x+1
    ud['iny'] = y+1

    ud['ymax'] = 8.0
    ud['rayleigh_forcing'] = False
    ud['mixed_run'] = True
    ud['do_advection'] = True
    ud['output_timesteps'] = False
    ud['stepmax'] = 10001


    for t in resol_t:
        ud['dtfixed0'] = t / t_ref
        ud['dtfixed']  = t / t_ref

        ud['rayleigh_forcing_fn'] = 'output_mark_wave_ensemble=1_%i_%i_36.000000_bottom_forcing_S%i.h5' %(x,int(4*y),t)
        ud['rayleigh_forcing_type'] = 'func'

    for om in omegas:
            ud['coriolis_strength'] = [0.0, 0.0, om]
            if om > 0:
                ud['aux'] = 'test_run_S%i_a05' %t
            else:
                if ud['mixed_run'] == True:
                    ud['aux'] = 'test_run_S%i_mix' %t
                else:
                    ud['aux'] = 'test_run_S%i_adv' %t

            print(ud)
            # run simulation
            rp.ud = json.dumps(ud)
            rp.dap = json.dumps(dap)
            rp.queue_run()

# path = '/home/ray/git-projects/RKLM_Reference/output_mark_wave/'
# # fn = 'output_mark_wave_ensemble=1_301_120_720.000000_rstrt_init_S600_a1.h5'
# fn = 'output_mark_wave_ensemble=1_301_120_720.000000_rstrt_init.h5'
# # name = '_002_after_full_step' # corresponding to output at 1800 s
# name = '_test' # corresponding to output at 1800 s

# # in dimensionless time
# ts = 18.0
# te = 720.0+6.0
# ti = 6.0

# rp.restart_set(path, fn, name, ts, te, ti)
# rp.restart_run()