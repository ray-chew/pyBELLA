# %%
from run import run_params as rp

from tests.diagnostics import compare_sol
import json

gen_targets    = False
updt_targets   = False

# define run instance
rp = rp()
rp.N = 1
ud = {}

if gen_targets:
    # set target user data parameters
    ud['output_type'] = 'target'
    ud['diag'] = False

    # generate target
    # define horizontal slice test case
    rp.tc = 'test_travelling_vortex'
    rp.ud = json.dumps(ud)
    rp.queue_run()

    # define vertical slice test case
    rp.tc = 'test_internal_long_wave'
    rp.ud = json.dumps(ud)
    rp.queue_run()

    # define test case
    rp.tc = 'test_lamb_wave'
    rp.ud = json.dumps(ud)
    rp.queue_run()

if updt_targets:
    diag = compare_sol("gen_target")
    diag.update_targets()




ud['output_type'] = 'test'
ud['diag'] = True

rp.tc = 'test_travelling_vortex'
rp.ud = json.dumps(ud)
rp.queue_run()


rp.tc = 'test_internal_long_wave'
rp.ud = json.dumps(ud)
rp.queue_run()


rp.tc = 'test_lamb_wave'
rp.ud = json.dumps(ud)
rp.queue_run()