from run import run_params as rp
import json

rp = rp()

# rp.N = 1
# rp.tc = 'tv'
# ud = {
#     'aux' : 'comp_bal',
#     # 'continuous_blending' : False
# }

# dap = {
#     'None' : None
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()

##########################################
#
# TV queue
#
##########################################

# rp.N = 10
# rp.tc = 'tv'
# ud = {
#     'aux' : 'debug_letkf_s10p',
#     'continuous_blending' : True
# }

# dap = {
#     # 'noise_percentage' : 0.1,
#     'None' : None
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()

# ud = {
#     'aux' : 'debug_letkf_s10p',
#     'continuous_blending' : False
# }

# rp.ud = json.dumps(ud)
# rp.queue_run()



##########################################
#
# SWE Queue, DA params
#
##########################################

rp.N = 10
rp.tc = 'swe_bal_vortex'


ud = {
    'aux' : 'debug_letkf_s05p',
    'continuous_blending' : True
}

dap = {
    'obs_attrs' : ['rhou', 'rhow'],
    'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
    'obs_frac' : 0.05
}

rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

rp.tc = 'swe_bal_vortex'

ud = {
    'aux' : 'debug_letkf_s05p',
    'continuous_blending' : False
}

dap = {
    'obs_attrs' : ['rhou', 'rhow'],
    'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
    'obs_frac' : 0.05
}

rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

ud = {
    'aux' : 'debug_letkf1_s025p',
    'continuous_blending' : True
}

dap = {
    'obs_attrs' : ['rhou', 'rhow'],
    'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
    'obs_frac' : 0.025
}

rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

ud = {
    'aux' : 'debug_letkf1_s25p',
    'continuous_blending' : False
}

dap = {
    'obs_attrs' : ['rhou', 'rhow'],
    'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
    'obs_frac' : 0.025
}

rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()


##########################################
#
# SWE Queue, Coriolis
#
##########################################

# rp.N = 1
# rp.tc = 'swe_bal_vortex'


# ud = {
#     'aux' : 'corr_1.0',
#     'continuous_blending' : False
# }

# dap = {
#     'None' : None,
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()

# rp.tc = 'swe_bal_vortex'

# ud = {
#     'aux' : 'debug_letkf_s50p',
#     'continuous_blending' : False
# }

# dap = {
#     'obs_attrs' : ['rhou', 'rhow'],
#     'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
#     'obs_frac' : 0.50
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()

# ud = {
#     'aux' : 'debug_letkf1_s10p',
#     'continuous_blending' : True
# }

# dap = {
#     'obs_attrs' : ['rhou', 'rhow'],
#     'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
#     'obs_frac' : 0.10
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()

# ud = {
#     'aux' : 'debug_letkf1_s10p',
#     'continuous_blending' : False
# }

# dap = {
#     'obs_attrs' : ['rhou', 'rhow'],
#     'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
#     'obs_frac' : 0.10
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()