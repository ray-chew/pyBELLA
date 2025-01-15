import logging
import copy

import numpy as np
import yaml

from ..utils.data_structures import DiagnosticState
from ..flow_solver.utils.solver_diagnostics import get_p_from_pressure_related_fields

from ..vis import (
    utils as vis_utils,
    plotting_tools as vis_pt
)


class CompareSol(object):
    def __init__(self, diag_state: DiagnosticState):
        self.diag_state = diag_state
        self.current_run = diag_state.test_name
        self.__init(diag_state)
        self.__get_tc()

    def update_targets(self):
        self.arr_dump = {}

        for tc_name, tc in self.tcs.items():
            tp = self.tps[tc_name]
            self.arr_dump[tp.name] = {}

            for attribute in tp.attributes:
                self.arr_dump[tp.name][attribute] = float(
                    self.__get_ens(tc, tp, attribute, summed=True)
                )

        with open("./src/tests/test_targets.yml", "a") as outfile:
            yaml.dump(self.arr_dump, outfile, default_flow_style=False)

    def test_do(self, mem, ud, plot=False):
        Sol = mem.sol
        mpv = mem.mpv

        self.__read_yaml()

        try:
            target_values = self.target[self.current_run]
        except:
            assert 0, "test %s has no target for comparison" % (self.current_run)

        if plot:
            self.__plot_comparison(mem, ud)

        for key, value in target_values.items():
            ref = value

            if key != "p2_nodes":
                test = getattr(Sol, key).astype("float32").sum()
            else:
                test = mpv.p2_nodes.astype("float32").sum()


            try:
                assert (
                    np.isclose(ref, test)
                ), "sum for attribute %s of %s changed with discrepancy:\n%.16f\n%.16f" % (
                    key,
                    self.current_run,
                    ref,
                    test,
                )
                logging.info(f"test passed for {key}")
            except AssertionError as e:
                logging.info(str(e))
                raise

        logging.info(f"""
        {'#' * 10}
        Test passed for {self.current_run}
        {'#' * 10}
        """.strip())

    def __init(self, ds: DiagnosticState
               ):

        tp = test_params(ds)

        self.tps = {
            ds.test_name: tp,
        }
        # self.tps = [tv_2D]

    def __get_tc(self):
        self.tcs = {}
        for test_name, test_param in self.tps.items():
            fn = test_param.fn + ".h5"
            tc = vis_utils.test_case(fn, test_param.dir, test_param.Nx, test_param.Ny, "")

            self.tcs[test_name] = tc

    def __read_yaml(self):
        with open("./src/tests/test_targets.yml", "r") as infile:
            self.target = yaml.safe_load(infile)

    def __plot_comparison(self, mem, ud):
        tc = self.tcs[self.current_run]
        tp = self.tps[self.current_run]

        ref_mem = copy.deepcopy(mem)
        for attribute in tp.attributes:

            if attribute != "p2_nodes":
                setattr(ref_mem.sol, attribute, self.__get_ens(tc, tp, attribute, summed=False))
            else:
                setattr(ref_mem.mpv, attribute, self.__get_ens(tc, tp, attribute, summed=False))

        for attribute in tp.attributes:
            arr_plots = []

            test_sol = self.__get_sol_for_comparison(mem, ud, attribute)
            ref_sol = self.__get_sol_for_comparison(ref_mem, ud, attribute)

            arr_plots.append([ref_sol, "ref"])
            arr_plots.append([test_sol, "test"])
            arr_plots.append([ref_sol - test_sol, "diff"])

            pl = vis_pt.plotter(arr_plots, ncols=3, figsize=(12, 3), sharey=False)
            _ = pl.plot(method="contour", lvls=None, suptitle=attribute)
            pl.img.savefig(tp.dir.replace("target", "test") + attribute + ".png")

        del ref_mem

    @staticmethod
    def __get_sol_for_comparison(mem, ud, attribute):
        Sol = mem.sol
        mpv = mem.mpv
        if attribute != "p2_nodes":
            test_sol = getattr(Sol, attribute).T
            if attribute != 'rho':
                rho = getattr(Sol, 'rho').T
                test_sol /= rho

                # if attribute == 'rhoY':
                #     test_sol -= mpv.HydroState.Y0[:,np.newaxis]

        else:
            test_sol = mpv.p2_nodes.T * ud.Msq
            test_sol -= mpv.HydroState_n.pi0[:,np.newaxis]
            # test_sol = get_p_from_pressure_related_fields(mem, ud).T

        return test_sol

    @staticmethod
    def __get_ens(tc, params, attribute, summed=True, normed=False):
        times = params.times
        l_typ = params.l_typ

        tags = tc.get_tag_dict()
        tag = "ic" if times[0] == 0.0 else tags[9]

        ens = tc.get_ensemble(
            times,
            1,
            attribute,
            "",
            label_type=l_typ,
            tag=tag,
            inner=False,
            get_fn=False,
            fn=params.fn + ".h5",
            load_ic=False,
            avg=True,
        )[0]

        if summed:
            return ens.sum()
        elif normed:
            return np.linalg.norm(ens)
        else:
            return ens


class test_params(object):
    def __init__(self, ds: DiagnosticState):

        self.name = ds.test_name
        self.dir = ds.path + ds.file_name + "/"
        self.fn = f"{ds.file_name}_{ds.Nx}_{ds.Ny}" 

        self.Nx = ds.Nx
        self.Ny = ds.Ny

        self.times = ds.steps
        self.l_typ = "WINDOW_STEP"

        self.attributes = ["rho", "rhou", "rhov", "rhow", "rhoY", "rhoX", "p2_nodes"]
