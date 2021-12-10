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
# 
#
##########################################
# 
gen_c6_swe_ic_half_full = True
gen_c6_swe_obs_truth = False
gen_c6_swe_enda = False

# Otherwise, if gen_all = True, generate all
# results
gen_all = False


# specify path to the observation directories 
# path_to_obs = '/srv/public/ray/'
path_to_obs = './'


if gen_c6_swe_ic_half_full or gen_all:
    ##########################################
    #
    #
    #
    ##########################################

    # No ensemble
    rp.N = 1 
    # SWE travelling vortex
    rp.tc = 'swe_bal_vortex'
    tout = [1.0]

    # JSON dumping does not accept ndarrays.
    # tout = tout.tolist()

    # lake run
    ud = {
        'aux' : 'lake_noib_imbal',
        # set pseudo-incompressible
        'is_compressible' : 0,
        'tout' : tout,
        'output_timesteps' : True
    }

    # No data assimilation.
    dap = {
        'None' : None,
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # simulation parameters for the 
    # SWE run without blending in the
    # initial time-step
    ud = {
        'aux' : 'imbal_noib',
        'tout' : tout,
        'output_timesteps' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()

    ##########################################

    # simulation parameters for the 
    # SWE run with blending in the
    # initial time-step

    # using pi-half
    ud = {
        'aux' : 'imbal_half',
        'initial_blending' : True,
        'tout' : tout,
        'output_timesteps' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()


    # using pi-full
    ud = {
        'aux' : 'imbal_full',
        'initial_blending' : True,
        'tout' : tout,
        'output_timesteps' : True,
        'blending_weight' : 1.0,
        'blending_type' : 'full'
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()


    ##########################################

    # simulation parameters for the 
    # SWE run without blending in the
    # initial time-step but with balanced IC
    ud = {
        'aux' : 'bal_noib',
        'tout' : tout,
        'output_timesteps' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    # rp.queue_run()



if gen_c6_swe_obs_truth or gen_all:
    ##########################################
    #
    # Generate observations and truths for the
    # shallow water vortices
    #
    ##########################################

    # No ensemble
    rp.N = 1 
    # Euler travelling vortex experiment
    rp.tc = 'swe_bal_vortex'

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



if gen_c6_swe_enda or gen_all:
    ##########################################
    #
    # Run simulations for the Euler vortex
    # ensembles. (section 6c)
    #
    ##########################################

    # Set ensemble with 10 members
    rp.N = 10
    # Euler travelling vortex experiment
    rp.tc = 'swe_bal_vortex'

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
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_obs.h5'
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # Generate ensemble with DA and without
    # blending
    ud = {
        'aux' : 'wda',
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
        'aux' : 'wda',
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


# if gen_5c1_euler_full or gen_all:
#     ##########################################
#     #
#     # Repeat the previous two experiments, but
#     # this time, assimilate all the quantities
#     #
#     ##########################################

#     # da_times = np.arange(5.0,10.5,0.5)/10.0
#     # da_times = np.array([0.25,10.0])
#     # da_times = np.arange(0.0,10.25,0.25)[1:]
#     da_times = np.arange(0.0,3.25,0.25)[1:]
#     da_times = np.around(da_times,3)
#     da_times = da_times.tolist()

#     # Set ensemble with 10 members
#     rp.N = 10
#     # Euler travelling vortex experiment
#     rp.tc = 'tv'
#     ud = {
#         'aux' : 'wda',
#         # Do blending for initial time-step
#         'initial_blending' : True
#     }

#     # Set the data assimilation parameters
#     dap = {
#         'da_times' : da_times,
#         # Assimilate all fields
#         'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
#         # Path to the generated observation
#         'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
#         'loc_setter' : (11,11)
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     # rp.queue_run()

#     ud = {
#         'aux' : 'wda',
#         # Do blending for initial time-step
#         'initial_blending' : True,
#         # Do blending after each assimilation
#         'continuous_blending' : True
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     rp.queue_run()

#     ##########################################



# if gen_loc_errors_tv or gen_all:
#     # Set ensemble with 10 members
#     rp.N = 10
#     # Euler travelling vortex experiment
#     rp.tc = 'tv'
#     ud = {
#         'aux' : 'wda_63',
#         # Do blending for initial time-step
#         'initial_blending' : True
#     }

#     # Set the data assimilation parameters
#     dap = {
#         # Assimilate all fields
#         'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
#         # Path to the generated observation
#         'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
#         'loc_setter' : (63,63)
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     # rp.queue_run()

#     ud = {
#         'aux' : 'wda_63',
#         # Do blending for initial time-step
#         'initial_blending' : True,
#         'continuous_blending' : True
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     # rp.queue_run()

    
#     ud = {
#         'aux' : 'wda_41',
#         # Do blending for initial time-step
#         'initial_blending' : True
#     }

#     # Set the data assimilation parameters
#     dap = {
#         # Assimilate all fields
#         'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
#         # Path to the generated observation
#         'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
#         'loc_setter' : (41,41)
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     rp.queue_run()

#     ud = {
#         'aux' : 'wda_41',
#         # Do blending for initial time-step
#         'initial_blending' : True,
#         'continuous_blending' : True
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     rp.queue_run()

    
#     ud = {
#         'aux' : 'wda_21',
#         # Do blending for initial time-step
#         'initial_blending' : True
#     }

#     # Set the data assimilation parameters
#     dap = {
#         # Assimilate all fields
#         'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
#         # Path to the generated observation
#         'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
#         'loc_setter' : (21,21)
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     rp.queue_run()

#     ud = {
#         'aux' : 'wda_21',
#         # Do blending for initial time-step
#         'initial_blending' : True,
#         'continuous_blending' : True
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     rp.queue_run()

#     ud = {
#         'aux' : 'wda_5',
#         # Do blending for initial time-step
#         'initial_blending' : True
#     }

#     # Set the data assimilation parameters
#     dap = {
#         # Assimilate all fields
#         'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
#         # Path to the generated observation
#         'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5',
#         'loc_setter' : (5,5)
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     rp.queue_run()

#     ud = {
#         'aux' : 'wda_5',
#         # Do blending for initial time-step
#         'initial_blending' : True,
#         'continuous_blending' : True
#     }

#     # run simulation
#     rp.ud = json.dumps(ud)
#     rp.dap = json.dumps(dap)
#     rp.queue_run()

if gen_c6_swe_ic_half_full == False and \
gen_c6_swe_obs_truth == False and \
gen_c6_swe_enda == False and \
gen_all == False:
    print("no results generated, all switches were set to False")
else:
    print("results_gen.py completed")
