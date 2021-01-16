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
gen_6b1 = True
gen_6b2 = False
# generate obs and truth for section 6c
gen_6c_obs_truth = False
# generate shallow water ensemble simulations
# for section 6c
gen_6c_swe = False
# generate Euler ensemble simulations for
# section 6c with assimilation of the
# momenta fields
gen_6c_euler_momenta = False
# generate Euler ensemble simulations for
# section 6c with assimilation of all
# fields
gen_6c_euler_full = False

# specify where the output directories are
path_to_obs = '/srv/public/ray/'
# path_to_obs = './'



if gen_6c_obs_truth:
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

    # run simulation
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

    ##########################################


    # shallow water vortex
    rp.tc = 'swe_bal_vortex' 

    # simulation parameters for the observation
    ud = {
        'aux' : 'obs',
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # simulation parameters for the truth
    ud = {
        'aux' : 'truth',
        # Do blending for initial time-step
        'initial_blending' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()



if gen_6c_swe:
    ##########################################
    #
    # Run simulations for the shallow water
    # vortex ensembles. (section 6c)
    #
    ##########################################

    # Set ensemble with 10 members
    rp.N = 10
    # Shallow water travelling vortex experiment
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
        'obs_path' : path_to_obs + 'output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_obs.h5',
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
        'obs_attrs' : ['rhou', 'rhow'],
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_obs.h5',
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # Generate ensemble with DA and with
    # blending
    ud = {
        'aux' : 'wdab',
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


if gen_6c_euler_momenta:
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
        'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5'
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
        'aux' : 'wdab',
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


if gen_6c_euler_full:
    ##########################################
    #
    # Repeat the previous two experiments, but
    # this time, assimilate all the quantities
    #
    ##########################################

    # Set ensemble with 10 members
    rp.N = 10
    # Euler travelling vortex experiment
    rp.tc = 'tv'
    ud = {
        'aux' : 'wda',
        # Do blending for initial time-step
        'initial_blending' : True
    }

    # Set the data assimilation parameters
    dap = {
        # Assimilate all fields
        'obs_attrs' : ['rho', 'rhou', 'rhov', 'rhoY', 'p2_nodes'],
        # Path to the generated observation
        'obs_path' : path_to_obs + 'output_travelling_vortex/output_travelling_vortex_ensemble=1_64_64_3.0_obs.h5'
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ud = {
        'aux' : 'wdab',
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




if gen_6b2:
    ##########################################
    #
    # Run simulations for the rising bubble
    # experiment (section 6b.2).
    #
    ##########################################

    # No ensemble
    rp.N = 1 
    # Rising bubble experiment
    rp.tc = 'rb'

    # simulation parameters for the pseudo-
    # incompressible rising bubble simulation
    ud = {
        'aux' : 'psinc',
        # set pseudo-incompressible simulation
        'is_compressible' : 0
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
    # compressible rising bubble without
    # blending in the initial time-step
    ud = {
        'aux' : 'comp_imbal',
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # simulation parameters for the 
    # compressible rising bubble with
    # blending in the initial time-step
    ud = {
        'aux' : 'comp_bal',
        'initial_blending' : True
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()




if gen_6b1:
    ##########################################
    #
    # Run simulations for the Euler and
    # shallow water vortices (section 6b.1).
    #
    ##########################################

    # No ensemble
    rp.N = 1 
    # Euler travelling vortex
    rp.tc = 'tv'
    # We run the simulation from t=0.0
    # to t=1.0, t non-dimensional.
    # Output at every 0.01 time units.
    tout = np.arange(0.0,1.01,0.01)[1:]

    # JSON dumping does not accept ndarrays.
    tout = tout.tolist()

    # simulation parameters for the pseudo-
    # incompressible run
    ud = {
        'aux' : 'psinc',
        # set pseudo-incompressible
        'is_compressible' : 0,
        'tout' : tout
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
    # compressible run without blending in the
    # initial time-step
    ud = {
        'aux' : 'comp_imbal',
        'tout' : tout
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # simulation parameters for the 
    # compressible run with blending in the
    # initial time-step
    ud = {
        'aux' : 'comp_imbal',
        'initial_blending' : True,
        'tout' : tout
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # Shallow water travelling vortex
    rp.tc = 'swe_bal_vortex'

    # simulation parameters for the pseudo-
    # incompressible run
    ud = {
        'aux' : 'psinc',
        # set pseudo-incompressible
        'is_compressible' : 0,
        'tout' : tout
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
    # compressible run without blending in the
    # initial time-step
    ud = {
        'aux' : 'comp_imbal',
        'tout' : tout
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()

    ##########################################

    # simulation parameters for the 
    # compressible run with blending in the
    # initial time-step
    ud = {
        'aux' : 'comp_imbal',
        'initial_blending' : True,
        'tout' : tout
    }

    # run simulation
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    rp.queue_run()


print("results_gen.py completed")