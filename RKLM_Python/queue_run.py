from run import run_params as rp
import json

rp = rp()

rp.N = 1
rp.tc = 'swe_dvortex'
ud = {
    # 'inx' : 64+1,
    # 'iny' : 128+1,
    'aux' : 'debug_incompressible'
}

dap = {
    'noise_percentage' : 0.1,

}
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

# rp.N = 1
# rp.tc = 'tv'
# ud = {
#     # 'inx' : 128+1,
#     # 'iny' : 128+1,
#     'aux' : 'debug_1',
#     'initial_blending' : False
# }

# dap = {
#     'noise_percentage' : 0.2,

# }
# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()
