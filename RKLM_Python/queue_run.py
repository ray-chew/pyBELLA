from run import run_params as rp
import json

rp = rp()

rp.N = 1
rp.tc = 'rb'
ud = {
    # 'inx' : 64+1,
    # 'iny' : 128+1,
    # 'aux' : 'comp_bal_noib'
    # 'aux' : 'comp_imbal_full'
    'aux' : 'debug_rkadv_comp_CFL_idt=test1'
    # 'None' : None
}

dap = {
    # 'noise_percentage' : 0.1,
    'None' : None
}

rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

# ud = {
#     'is_compressible' : 1,
#     'compressibility' : 1.0,
#     'aux' : 'debug_rkadv_comp_CFL_idt=1'
# }

# rp.ud = json.dumps(ud)
# rp.queue_run()

# ud = {
#     'is_compressible' : 1,
#     'compressibility' : 1.0,
#     'initial_blending' : True,
#     'aux' : 'debug_rkadv_comp_CFL_idt=1'
# }

# rp.ud = json.dumps(ud)
# rp.queue_run()
