from run import run_params as rp
import json
import numpy as np

rp = rp()
path_to_obs = './'

rp.N = 1
rp.tc = 'swe_bal_vortex'

dap = {
    'None': None
}

rp.dap = json.dumps(dap)

### Travelling SWE Vortex
tout = np.arange(0.0, 1.01, 0.01)[1:]
tout = tout.tolist()


## 64x64
ud = {
    'aux' : 'debug_dsi_tra',
    'inx' : 64+1,
    'inz' : 64+1,
    'tout' : tout
}

# run simulation
rp.ud = json.dumps(ud)
rp.queue_run()


## 128x128
ud = {
    'aux' : 'debug_dsi_tra',
    'inx' : 128+1,
    'inz' : 128+1,
    'tout' : tout
}

# run simulation
rp.ud = json.dumps(ud)
rp.queue_run()


## 256x256
ud = {
    'aux' : 'debug_dsi_tra',
    'inx' : 256+1,
    'inz' : 256+1,
    'tout' : tout
}

# run simulation
rp.ud = json.dumps(ud)
rp.queue_run()



### Stationary SWE Vortex
## 64x64
ud = {
    'aux' : 'debug_dsi_sta',
    'inx' : 64+1,
    'inz' : 64+1,
    'tout' : tout,
    'u_wind_speed' : 0.0,
    'w_wind_speed' : 0.0
}

# run simulation
rp.ud = json.dumps(ud)
rp.queue_run()


## 128x128
ud = {
    'aux' : 'debug_dsi_sta',
    'inx' : 128+1,
    'inz' : 128+1,
    'tout' : tout,
    'u_wind_speed' : 0.0,
    'w_wind_speed' : 0.0
}

# run simulation
rp.ud = json.dumps(ud)
rp.queue_run()


## 256x256
ud = {
    'aux' : 'debug_dsi_sta',
    'inx' : 256+1,
    'inz' : 256+1,
    'tout' : tout,
    'u_wind_speed' : 0.0,
    'w_wind_speed' : 0.0
}

# run simulation
rp.ud = json.dumps(ud)
rp.queue_run()





