# import numpy as np
import yaml

# import sys
# sys.path.append('..')

from visualiser_debugger import utils
# import plotting_tools as pt


class compare_sol(object):
    def __init__(self, current_run):
        self.current_run = current_run

    def update_target(self):
        self.__init()
        self.__get_tc()

        self.arr_dump = {}

        for idx, tc in enumerate(self.tcs):
            tp = self.tps[idx]
            self.arr_dump[tp.name] =  {}

            for attribute in tp.attributes:
                self.arr_dump[tp.name][attribute] =                float(self.__get_ens(tc, tp, attribute))

        with open('targets.yml', 'w') as outfile:
            yaml.dump(self.arr_dump, outfile, default_flow_style=False)


    def compare(self, Sol, p2n):
        self.__read_yaml()

        try: 
            target_values = self.target[self.current_run]
        except:
            assert 0, "test %s has no target for comparison" %(self.current_run)


        for key, value in vars(Sol).items():
            assert target_values[key] == value, "test failed for attribute %s of %s" %(key, self.current_run)

        assert target_values['p2_nodes'] == p2n, "test failed for attribute p2_nodes of %s" %(self.current_run)


    def __read_yaml(self):
        with open('targets.yml', 'r') as infile:
            self.target = yaml.safe_load(infile)


    def __init(self):
        path = './RKLM_Python/tests/targets/' 

        tv_2D = test_params('test_travelling_vortex', path, 'target_travelling_vortex', 64, 64, [100])
        igw   = test_params('test_internal_long_wave', path, 'target_internal_long_wave', 301, 10, [30])
        lmbw  = test_params('test_lamb_wave', path, 'target_lamb_wave', 301, 120, [30])

        # self.tps = [tv_2D, igw, lmbw]
        self.tps = [tv_2D]


    def __get_tc(self):
        self.tcs = []
        for test_param in self.tps:
            tc = utils.test_case(test_param.fn, test_param.dir , test_param.Nx, test_param.Ny, '')

            self.tcs.append(tc)


    @staticmethod
    def __get_ens(tc, params, attribute):
        times = params.times
        l_typ = params.l_typ

        tags = tc.get_tag_dict()
        tag = 'ic' if times[0] == 0.0 else tags[9]

        ens = tc.get_ensemble(times, 1, attribute, '', label_type=l_typ, tag=tag, inner=True, get_fn=False, fn=params.fn, load_ic=False)

        return ens.sum()


class test_params(object):
    def __init__(self, name, path, fn, Nx, Ny, times):
        
        self.name = name
        self.dir = path
        self.fn = '%s_%i_%i.h5' %(fn, Nx, Ny)

        self.Nx = Nx
        self.Ny = Ny

        self.times = times
        self.l_typ = 'WINDOW_STEP'

        self.attributes = [
            'rho',
            'rhou',
            'rhov',
            'rhoY',
            'p2_nodes'
        ]

