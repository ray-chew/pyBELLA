from run import run_params as rp
import json
import numpy as np

rp = rp()
path_to_obs = './'

rp.N = 1

rp.tc = 'mark'

t_ref = 100.0
omega = 7.292 * 1e-5
resol_x = [151,301,601]
resol_y = [15,30,60]
# resol_x = [301]
# resol_y = [30]
resol_t = [200,400,600,800,1000,1200,1400,1600]
# resol_t = [200,400,600,800]
# resol_t = [400,600]
omegas = [0.0, 2.0 * omega * t_ref]
# omegas = [2.0 * omega * t_ref]


ud = {}
dap = {
    'None' : None,
}

for x,y in zip(resol_x,resol_y):
    ud['inx'] = x+1
    ud['iny'] = y+1

    for t in resol_t:
        ud['dtfixed0'] = t / t_ref
        ud['dtfixed']  = t / t_ref

        for om in omegas:
            ud['coriolis_strength'] = [0.0, 0.0, om]
            if om > 0:
                ud['aux'] = 'bdl_run_S%i_a05' %t
            else:
                ud['aux'] = 'bdl_run_S%i_noom_a05' %t

            print(ud)
            # run simulation
            rp.ud = json.dumps(ud)
            rp.dap = json.dumps(dap)
            rp.queue_run()