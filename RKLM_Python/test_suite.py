# %%
from run import run_params as rp

from tests.diagnostics import compare_sol
import json

%load_ext autoreload

import os

# os.chdir(os.getcwd())
print(os.getcwd())


# define run instance
rp = rp()
rp.N = 1

ud = {}

# # define horizontal slice test case
ud['output_type'] = 'target'
rp.tc = 'test_travelling_vortex'
rp.ud = json.dumps(ud)
rp.queue_run()

# # define vertical slice test case
# rp.tc = 'test_intertia_gravity_wave'
# rp.queue_run()

# # define test case
# rp.tc = 'test_lamb_coriolis'
# rp.queue_run()


# %%
%autoreload
diag = compare_sol("test_travelling_vortex")

diag.update_target()

# %%
