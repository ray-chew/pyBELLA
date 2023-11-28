# %%
from run import run_params as rp

from tests.diagnostics import compare_sol
import json

# %load_ext autoreload
# %autoreload
# %matplotlib inline
import os

# os.chdir(os.getcwd())
print(os.getcwd())


gen_targets = True

# define run instance
rp = rp()
rp.N = 1
ud = {}

if gen_targets:
    # generate target
    # define horizontal slice test case
    rp.tc = 'test_travelling_vortex'
    ud['output_type'] = 'target'
    ud['diag'] = False
    rp.ud = json.dumps(ud)
    rp.queue_run()

    diag = compare_sol("test_travelling_vortex")
    diag.update_target()
    
    # # define vertical slice test case
    # rp.tc = 'test_intertia_gravity_wave'
    # rp.queue_run()

    # # define test case
    # rp.tc = 'test_lamb_coriolis'
    # rp.queue_run()

rp.tc = 'test_travelling_vortex'
ud['output_type'] = 'test'
ud['diag'] = True
rp.ud = json.dumps(ud)
rp.queue_run()