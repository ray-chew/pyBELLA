import numpy as np
import yaml

import sys
sys.path.append('..')

from visualiser_debugger import utils
# import plotting_tools as pt


class diagnostics(object):
    def __init__(self, current_run):
        self.current_run = current_run


    def __init(self):
        path = './targets/' 

        tv_2D = test_params(path, 'target_travelling_vortex', 64, 64, [100])
        igw   = test_params(path, 'target_internal_long_wave', 301, 10, [30])
        lmbw  = test_params(path, 'target_lamb_wave', 301, 120, [30])

        self.test_params = [tv_2D, igw, lmbw]


    def update_target(self):
        with open('data.yml', 'w') as outfile:
            yaml.dump(data, outfile, default_flow_style=False)

    


    def __get_tc(self):
        self.tcs = []
        for test_param in self.test_params:
            tc = utils.test_case(test_param.fn, test_param.dir , test_param.Nx, test_param.Ny, '')

            self.tcs.append(tc)



    def get_ens(tc, params, attribute):
        tags = tc.get_tag_dict()
        tag = 'ic' if times[0] == 0.0 else tags[9]

        times = params.times
        l_typ = params.l_typ

        diff = False

        ens = tc.get_ensemble(times, 1, attribute, '', label_type=l_typ, avg=True, diff=diff, tag=tag, inner=True)[1]
        # rho = tc.get_ensemble(times, 1, 'rho', '', label_type=l_typ, avg=True, diff=diff, tag=tag, inner=True)[1]

        # ens = ens.T / rho.T
        return attribute, ens


class test_params(object):
    def __init__(self, path, fn, Nx, Ny, times):
        
        self.dir = path
        self.fn = fn

        self.Nx = Nx
        self.Ny = Ny

        self.times = times
        self.l_typ = 'WINDOW_STEP'

        self.attributes = [
            'rho',
            'rhou',
            'rhov',
            'rhoY',
            'pi'
        ]

