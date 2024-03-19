# %%
from run import run_params as rp
from tests.diagnostics import compare_sol
import json

# The first run MUST have gen_targets = True
# This is to generate the target outputs
# so that a comparison can be made
gen_targets = False

# Only set this to True if the diagnostic targets
# are to be updated.
updt_targets = False

# define run instance
rp = rp()
rp.N = 1
ud = {}

if gen_targets:
    # set target user data parameters
    ud["output_type"] = "target"
    ud["diag"] = False

    # generate target
    # define horizontal slice test case
    rp.tc = "test_travelling_vortex"
    rp.ud = json.dumps(ud)
    rp.queue_run()

    # define vertical slice test case
    rp.tc = "test_internal_long_wave"
    rp.ud = json.dumps(ud)
    rp.queue_run()

    # define test case
    rp.tc = "test_lamb_wave"
    rp.ud = json.dumps(ud)
    rp.queue_run()

if updt_targets:
    diag = compare_sol("gen_target")
    diag.update_targets()


ud["output_type"] = "test"
# Do diagnostics
ud["diag"] = True

rp.tc = "test_travelling_vortex"
rp.ud = json.dumps(ud)
rp.queue_run()

rp.tc = "test_internal_long_wave"
rp.ud = json.dumps(ud)
rp.queue_run()

rp.tc = "test_lamb_wave"
rp.ud = json.dumps(ud)
rp.queue_run()
