from run import run_params as rp
import json
import numpy as np

rp = rp()
path_to_obs = './'

##########################################
#
# Rising bubble ensemble generation
#
##########################################

# Generate obs and truth (which are the same here)
# rp.N = 1
rp.tc = 'rb'

ud = {
    'aux' : 'truth_CFLfixed',
    'initial_blending' : True
}

dap = {
    'None' : None,
}

# run simulation
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

##########################################

# Generate ensemble simulations
rp.N = 10

ud = {
    'aux' : 'noda_CFLfixed',
    'initial_blending' : True
}

dap = {
    'da_times' : [],
    'obs_path' : path_to_obs + 'output_rising_bubble/output_rising_bubble_ensemble=1_160_80_10.0_truth_CFLfixed_ib-0.h5'
}

# run simulation
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

##########################################

# Generate ensemble with DA and without
# blending
ud = {
    'aux' : 'wda_CFLfixed',
    # Do blending for initial time-step
    'initial_blending' : True
}

da_times = np.arange(5.0,10.5,0.5)
da_times = da_times.tolist()

# Set the data assimilation parameters
dap = {
    'da_times' : da_times,
    # Assimilate the momentum fields
    'obs_attrs' : ['rhou', 'rhov'],
    # Path to the generated observation
    'obs_path' : path_to_obs + 'output_rising_bubble/output_rising_bubble_ensemble=1_160_80_10.0_truth_CFLfixed_ib-0.h5'
}

# run simulation
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

##########################################

# Generate ensemble with DA and with
# blending
ud = {
    'aux' : 'wda_CFLfixed',
    # Do blending for initial time-step
    'initial_blending' : True,
    # Do blending after each assimilation
    'continuous_blending' : True
}

# We do not have to change the DA parameters
# so we use 'dap' dictionary from before.

# run simulation
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()