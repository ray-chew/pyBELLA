from run import run_params as rp
import json

rp = rp()

##########################################
#
# TV queue, obs RMSEs
#
##########################################

rp.N = 10
rp.tc = 'tv_neg'
ud = {
    'aux' : 'debug_neg_noda',
    'continuous_blending' : False,
    'initial_blending' : True
}

dap = {
    'da_times' : [],
}
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

ud = {
    'aux' : '',
    'continuous_blending' : False,
    'initial_blending' : True
}

dap = {
    'none' : None,
}
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

ud = {
    'aux' : '',
    'continuous_blending' : True,
    'initial_blending' : True
}

dap = {
    'none' : None,
}

rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()


##########################################
#
# TV queue, obs RMSEs
#
##########################################

# rp.N = 10
# rp.tc = 'tv'
# ud = {
#     'aux' : 'debug_pos_s10p_n15p',
#     'continuous_blending' : False
# }

# dap = {
#     'noise_percentage' : 0.15,
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()

# ud = {
#     'aux' : 'debug_pos_s10p_n20p',
#     'continuous_blending' : False
# }

# dap = {
#     'noise_percentage' : 0.20,
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()


# ud = {
#     'aux' : 'debug_pos_s10p_n30p',
#     'continuous_blending' : False
# }

# dap = {
#     'noise_percentage' : 0.30,
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()


# ud = {
#     'aux' : 'debug_pos_s10p_n35p',
#     'continuous_blending' : False
# }

# dap = {
#     'noise_percentage' : 0.35,
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()


# ud = {
#     'aux' : 'debug_pos_s10p_n40p',
#     'continuous_blending' : False
# }

# dap = {
#     'noise_percentage' : 0.40,
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()


# ud = {
#     'aux' : 'debug_pos_s10p_n45p',
#     'continuous_blending' : False
# }

# dap = {
#     'noise_percentage' : 0.45,
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()



##########################################
#
# SWE Queue, DA params
#
##########################################

# rp.N = 10
# rp.tc = 'swe_bal_vortex'


# ud = {
#     'aux' : 'debug_letkf_s05p',
#     'continuous_blending' : True
# }

# dap = {
#     'obs_attrs' : ['rhou', 'rhow'],
#     'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
#     'obs_frac' : 0.05
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()

# rp.tc = 'swe_bal_vortex'

# ud = {
#     'aux' : 'debug_letkf_s05p',
#     'continuous_blending' : False
# }

# dap = {
#     'obs_attrs' : ['rhou', 'rhow'],
#     'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
#     'obs_frac' : 0.05
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()

# ud = {
#     'aux' : 'debug_letkf1_s025p',
#     'continuous_blending' : True
# }

# dap = {
#     'obs_attrs' : ['rhou', 'rhow'],
#     'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
#     'obs_frac' : 0.025
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()

# ud = {
#     'aux' : 'debug_letkf1_s25p',
#     'continuous_blending' : False
# }

# dap = {
#     'obs_attrs' : ['rhou', 'rhow'],
#     'obs_path' : './output_swe_vortex/output_swe_vortex_ensemble=1_64_1_64_3.0_comp_1.0_pp_tra_truth_ip.h5',
#     'obs_frac' : 0.025
# }

# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()


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

# rp.N = 1
# rp.tc = 'swe_bal_vortex'
# ud = {
#     'inx' : 128+1,
#     'inz' : 128+1,
#     'aux' : 'debug_dsi_da'
# }
# 
# dap = {
#     'noise_percentage' : 0.1,
# 
# }
# rp.ud = json.dumps(ud)
# rp.dap = json.dumps(dap)
# rp.queue_run()
# 
# rp.N = 1
# rp.tc = 'swe_bal_vortex'
# ud = {
#     'inx' : 256+1,
#     'inz' : 256+1,
#     'aux' : 'debug_dsi_da'
# }
# 
# dap = {
#     'noise_percentage' : 0.1,
# 
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