from run import run_params as rp
import json

rp = rp()


rp.N = 1
rp.tc = 'tv'
ud = {
    'inx' : 128+1,
    'iny' : 128+1,
    'aux' : 'debug'
}
rp.ud = json.dumps(ud)

rp.queue_run()

