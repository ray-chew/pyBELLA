from run import run_params as rp
import json

# define run instance
rp = rp()
rp.N = 1

# # define horizontal slice test case
rp.tc = 'test_travelling_vortex'
rp.queue_run()

# # define vertical slice test case
rp.tc = 'test_intertia_gravity_wave'
rp.queue_run()

# define test case
rp.tc = 'test_lamb_coriolis'
rp.queue_run()

