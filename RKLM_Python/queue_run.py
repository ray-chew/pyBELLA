from run import run_params as rp
import json

rp = rp()

rp.N = 1
rp.tc = 'swe_dvortex'
ud = {
    'inx' : 64+1,
    'inz' : 64+1,
    'aux' : 'debug_dsi_da_incompressible'
}

dap = {
    'noise_percentage' : 0.1,

}

rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()
rp.N = 1
rp.tc = 'swe_dvortex'
ud = {
    'inx' : 128+1,
    'inz' : 128+1,
    'aux' : 'debug_dsi_da_incompressible'
}

dap = {
    'noise_percentage' : 0.1,

}
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

rp.N = 1
rp.tc = 'swe_dvortex'
ud = {
    'inx' : 256+1,
    'inz' : 256+1,
    'aux' : 'debug_dsi_da_incompressible'
}

dap = {
    'noise_percentage' : 0.1,

}
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()
