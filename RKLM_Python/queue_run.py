from run import run_params as rp
import json

rp = rp()

rp.N = 3
rp.tc = 'tv'
ud = {
    # 'inx' : 128+1,
    # 'iny' : 128+1,
    'aux' : 'debug'
}

dap = {
    'noise_percentage' : 0.1,

}
rp.ud = json.dumps(ud)
rp.dap = json.dumps(dap)
rp.queue_run()

