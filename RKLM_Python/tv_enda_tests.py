##########################################
#
# This script generates the results from
# the RKLM-Py code for the publication
# "Balanced Local Data Assimilation with
# a Blended Numerical Model for Horizontal
# Dynamics".
# 
# Refer to the quickstart document for
# more details on how to get the code
# running.
#
##########################################

##########################################
#
# A note on reproducibility:
# The initial center position of the
# vortex as well as the ensemble spread
# are generated randomly. These random
# functions are seeded for reproducibility
# . Nevertheless, seeded random functions
# may produce slightly different results
# across platforms and architectures.
#
##########################################

from run import run_params as rp
import json
import numpy as np

rp = rp()

##########################################
#
# This script is divided into three parts.
# Each part generates results for the
# sections 6b.1), 6b.2) and 6c.
# Set the section to True to generate the
# results.
#
##########################################
# generate obs and truth for section 5b
gen_5b_obs_truth_euler = False
# generate Euler ensemble simulations for
# section 5c1 with assimilation of all
# fields
generate_noda = False
gen_5c1_euler_full = False
gen_5c1_euler_momenta = False
gen_loc_errors_tv = True

# Otherwise, if gen_all = True, generate all
# results
gen_all = False

# specify path to the observation directories 
# path_to_obs = '/srv/public/ray/'
path_to_obs = './'

enda_sfx = 'wda_varcov'

if gen_5b_obs_truth_euler or gen_all:
    ##########################################
    #
    # Generate observations and truths for the
    # Euler and shallow water vortices
    # (section 6c)
    #
    ##########################################

    # No ensemble
    rp.N = 1 
    # Euler travelling vortex experiment
    rp.tc = 'tv'

    # simulation parameters for the observation
    ud = {
        'aux' : 'obs',
    }

    # data assimilation parameters for the
    # observation. No DA.
    dap = {
        'None' : None,
    }

    # # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # simulation parameters for the truth
    ud = {
        'aux' : 'truth',
        # Do blending for the initial time-step
        'initial_blending' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()


if gen_5c1_euler_full or gen_all:
    ##########################################
    #
    # Repeat the previous two experiments, but
    # this time, assimilate all the quantities
    #
    ##########################################

    da_times = np.arange(0.0,10.25,0.25)[1:]
    da_times = np.arange(0.0,3.25,0.25)[1:]
    da_times = np.around(da_times,3)
    da_times = da_times.tolist()

    # Set ensemble with 10 members
    rp.N = 10
    # Euler travelling vortex experiment
    rp.tc = 'tv'

    # Generate ensemble with no DA
    ud = {
        'aux' : 'noda',
        # Do blending for initial time-step
        'initial_blending' : True
    }

    # For the data assimilation parameters,
    # do not do DA at any time-point.
    dap = {
        'da_times' : [],
        'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
        'loc_setter' : (11, 11),
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5'
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    if generate_noda:
        rp.queue_run()

    ud = {
        'aux' : '%s' %enda_sfx,
        # Do blending for initial time-step
        'initial_blending' : True
    }

    # Set the data assimilation parameters
    u_ref = 10000.0 / 100.0
    T_ref = 300.0
    p_ref = 1e+5
    T = 1.0 / T_ref
    p = (50 / p_ref)
    rho = p / T
    pi = p**0.4
    v = 4.0 / u_ref

    dap = {
        'da_times' : da_times,
        # Assimilate all fields
        'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
        'loc_setter' : (11, 11),
        # using variance: 1K for temperature, 4ms^-1 for velocity, 50Pa for pressure
        # 'sd_setter' : [rho**0.5, (rho*v)**0.5, (rho*v)**0.5, p**0.5, pi**0.5],
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ud = {
        'aux' : '%s' %enda_sfx,
        # Do blending for initial time-step
        'initial_blending' : True,
        # Do blending after each assimilation
        'continuous_blending' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################


if gen_5c1_euler_momenta or gen_all:
    ##########################################
    #
    # Run simulations for the Euler vortex
    # ensembles. (section 6c)
    #
    ##########################################

    # Set ensemble with 10 members
    rp.N = 10
    # Euler travelling vortex experiment
    rp.tc = 'tv'

    ##########################################

    # Generate ensemble with DA and without
    # blending
    ud = {
        'aux' : '%s' %enda_sfx,
        # Do blending for initial time-step
        'initial_blending' : True
    }

    # Set the data assimilation parameters
    dap = {
        # Assimilate the momentum fields
        'obs_attrs' : ['rhou', 'rhov'],
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5'
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # Generate ensemble with DA and with
    # blending
    ud = {
        'aux' : '%s' %enda_sfx,
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




if gen_loc_errors_tv or gen_all:
    da_times = np.arange(0.0,3.25,0.25)[1:]
    da_times = np.around(da_times,3)
    da_times = da_times.tolist()

    # Set ensemble with 10 members
    rp.N = 10
    # Euler travelling vortex experiment
    rp.tc = 'tv'
    ud = {
        'aux' : 'wda_63',
        # Do blending for initial time-step
        'initial_blending' : True
    }

    # Set the data assimilation parameters
    dap = {
        'da_times' : da_times,
        # Assimilate all fields
        'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
        'loc_setter' : (63,63)
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()

    ud = {
        'aux' : 'wda_63',
        # Do blending for initial time-step
        'initial_blending' : True,
        'continuous_blending' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    
    ud = {
        'aux' : 'wda_31',
        # Do blending for initial time-step
        'initial_blending' : True
    }

    # Set the data assimilation parameters
    dap = {
        'da_times' : da_times,
        # Assimilate all fields
        'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
        'loc_setter' : (31,31)
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()

    ud = {
        'aux' : 'wda_31',
        # Do blending for initial time-step
        'initial_blending' : True,
        'continuous_blending' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()

    
    ud = {
        'aux' : 'wda_15',
        # Do blending for initial time-step
        'initial_blending' : True
    }

    # Set the data assimilation parameters
    dap = {
        'da_times' : da_times,
        # Assimilate all fields
        'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
        'loc_setter' : (15,15)
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()

    ud = {
        'aux' : 'wda_15',
        # Do blending for initial time-step
        'initial_blending' : True,
        'continuous_blending' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()

    ud = {
        'aux' : 'wda_7',
        # Do blending for initial time-step
        'initial_blending' : True
    }

    # Set the data assimilation parameters
    dap = {
        'da_times' : da_times,
        # Assimilate all fields
        'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
        'loc_setter' : (7,7)
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()

    ud = {
        'aux' : 'wda_7',
        # Do blending for initial time-step
        'initial_blending' : True,
        'continuous_blending' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()



if gen_5b_obs_truth_euler == False and \
gen_5c1_euler_full == False and \
gen_5c1_euler_momenta == False and \
gen_loc_errors_tv == False and \
gen_all == False:
    print("no results generated, all switches were set to False")
else:
    print("tv_enda_tests.py completed")