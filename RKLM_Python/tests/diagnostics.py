# import numpy as np
import yaml

import sys
sys.path.append("..") 
sys.path.append(".")


from visualiser_debugger import utils
from visualiser_debugger import plotting_tools as pt
# import plotting_tools as pt

class compare_sol(object):
    def __init__(self, current_run):
        self.current_run = current_run
        self.__init()
        self.__get_tc()

    def update_target(self):
        self.arr_dump = {}

        for idx, tc in enumerate(self.tcs):
            tp = self.tps[idx]
            self.arr_dump[tp.name] =  {}

            for attribute in tp.attributes:
                self.arr_dump[tp.name][attribute] =                float(self.__get_ens(tc, tp, attribute, summed=True))

        with open('targets.yml', 'w') as outfile:
            yaml.dump(self.arr_dump, outfile, default_flow_style=False)


    def do(self, Sol, p2n, plot=False):
        self.__read_yaml()

        try: 
            target_values = self.target[self.current_run]
        except:
            assert 0, "test %s has no target for comparison" %(self.current_run)

        if plot:
            self.__plot_comparison(Sol, p2n)

        for key, value in target_values.items():
            ref = value

            if key != 'p2_nodes':
                test = getattr(Sol, key).astype('float32').sum()
            else:
                test = p2n.astype('float32').sum()

            assert  ref == test , "test failed for attribute %s of %s with discrepancy:\n%.6f\n%.6f" %(key, self.current_run, ref, test)


    def __init(self):
        path = './' 

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


    def __read_yaml(self):
        with open('targets.yml', 'r') as infile:
            self.target = yaml.safe_load(infile)


    def __plot_comparison(self, Sol, p2n):
        for idx, tc in enumerate(self.tcs):
            tp = self.tps[idx]

            for attribute in tp.attributes:
                arr_plots = []

                ref_sol = self.__get_ens(tc, tp, attribute)

                if attribute != 'p2_nodes':
                    test_sol = getattr(Sol, attribute)
                else:
                    test_sol = p2n

                arr_plots.append([ref_sol, "ref"])
                arr_plots.append([test_sol, "test"])
                arr_plots.append([ref_sol - test_sol, "diff"])

                pl = pt.plotter(arr_plots, ncols=3, figsize=(12,3),sharey=False)
                _ = pl.plot(method='contour', lvls=None, suptitle=attribute)
                pl.img.show()

    @staticmethod
    def __get_ens(tc, params, attribute, summed=False):
        times = params.times
        l_typ = params.l_typ

        tags = tc.get_tag_dict()
        tag = 'ic' if times[0] == 0.0 else tags[9]

        ens = tc.get_ensemble(times, 1, attribute, '', label_type=l_typ, tag=tag, inner=False, get_fn=False, fn=params.fn, load_ic=False, avg=True)[0]

        if summed:
            return ens.sum()
        else:
            return ens


class test_params(object):
    def __init__(self, name, path, fn, Nx, Ny, times):
        
        self.name = name
        self.dir = path + fn + '/'
        self.fn = '%s_%i_%i.h5' %(fn, Nx, Ny)

        self.Nx = Nx
        self.Ny = Ny

        self.times = times
        self.l_typ = 'WINDOW_STEP'

        self.attributes = [
            'rho',
            'rhou',
            'rhov',
            'rhow',
            'rhoY',
            'rhoX',
            'p2_nodes'
        ]

