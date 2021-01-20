from run import run_params as rp
import json
import numpy as np

rp = rp()
path_to_obs = './'

rp.N = 1

# Euler travelling vortex experiment
rp.tc = 'tv'

# # simulation parameters for the observation
ud = {
    'aux' : 'obs_homen',
}

# # data assimilation parameters for the
# # observation. No DA.
dap = {
    'None' : None,
}

# # run simulation
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

# ##########################################

# # simulation parameters for the truth
ud = {
    'aux' : 'truth_homen',
    # Do blending for the initial time-step
    'initial_blending' : True
}

# run simulation
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

##########################################

##########################################
#
# Run simulations for the Euler vortex
# ensembles.
#
##########################################

# Set ensemble with 10 members
rp.N = 10
# Euler travelling vortex experiment
rp.tc = 'tv'

# Generate ensemble with no DA
ud = {
    'aux' : 'noda_homen',
    # Do blending for initial time-step
    'initial_blending' : True
}

# For the data assimilation parameters,
# do not do DA at any time-point.
dap = {
    'da_times' : [],
    # Path to the generated observation
    'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs_homen.h5'
}

# run simulation
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

##########################################

# Generate ensemble with DA and without
# blending
ud = {
    'aux' : 'wda_homen',
    # Do blending for initial time-step
    'initial_blending' : True
}

# Set the data assimilation parameters
dap = {
    # Assimilate the momentum fields
    'obs_attrs' : ['rhou', 'rhov'],
    # Path to the generated observation
    'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs_homen.h5'
}

# run simulation
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

##########################################

# Generate ensemble with DA and with
# blending
ud = {
    'aux' : 'wda_homen',
    # Do blending for initial time-step
    'initial_blending' : True,
    # Do blending after each assimilation
    'continuous_blending' : True
}

# We do not have to set the DA parameters
# so we use 'dap' dictionary from before.

# run simulation
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()





