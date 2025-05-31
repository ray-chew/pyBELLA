import logging
import copy

import numpy as np
import yaml

from ..utils.data_structures import DiagnosticState
from ..flow_solver.utils.solver_diagnostics import get_p_from_pressure_related_fields

from ..vis import utils as vis_utils, plotting_tools as vis_pt


class CompareSol(object):
    def __init__(self, diag_state: DiagnosticState):
        self.diag_state = diag_state
        self.current_run = diag_state.test_name
        self.plot = diag_state.plot_compare
        self.tolerances = diag_state.tolerances
        self.time_increment = diag_state.time_increment
        self.__init(diag_state)
        self.__get_tc()

    def update_targets(self):
        self.arr_dump = {}

        for tc_name, tc in self.tcs.items():
            tp = self.tps[tc_name]
            dump_name = tp.name.replace("target","test")
            self.arr_dump[dump_name] = {}

            for attribute in tp.attributes:
                arr = self.__get_ens(tc, tp, attribute, time_increment=self.time_increment, summed=False)
                self.arr_dump[dump_name][attribute] = float(
                    arr.sum()
                )

                if self.plot:
                    # vis_pt.plotter accepts a list of tuples with plot and panel title.
                    pl = vis_pt.plotter([(arr.T, "ref"), ], ncols=1, figsize=(4, 3), sharey=False)
                    _ = pl.plot(method="contour", lvls=None, suptitle=attribute)
                    pl.img.savefig(tp.dir + attribute + ".png")

        with open("./src/tests/test_targets.yml", "a") as outfile:
            yaml.dump(self.arr_dump, outfile, default_flow_style=False)

    def test_do(self, mem, ud):
        tc = self.tcs[self.current_run]
        tp = self.tps[self.current_run]

        # populate reference arrays
        ref_mem = copy.deepcopy(mem)
        for attribute in tp.attributes:
            try:
                ref_data = self.__get_ens(tc, tp, attribute, time_increment=False, summed=False)
                if attribute != "p2_nodes":
                    setattr(ref_mem.sol, attribute, ref_data)
                else:
                    if self.time_increment:
                        # in the case where we are interested in the time increment
                        # of the pressure-related fields, we need to load the previous
                        # time step from the output file.
                        tc_test = copy.deepcopy(tc)
                        tc_test.base_fn = tc_test.base_fn.replace("target", "test")
                        tc_test.py_dir = tc_test.py_dir.replace("target", "test")
                        data = self.__get_ens(tc_test, tp, attribute, time_increment=self.time_increment, summed=False)
                        setattr(mem.mpv, attribute, data)
                        ref_data = self.__get_ens(tc, tp, attribute, time_increment=self.time_increment, summed=False)
                    setattr(ref_mem.mpv, attribute, ref_data)
            except Exception as e:
                raise AssertionError(f"test {self.current_run} has no target for comparison: {e}")

        if self.plot:
            self.__plot_comparison(mem, ref_mem, ud)

        for attribute in tp.attributes:
            test = self.__get_sol_for_comparison(mem, ud, attribute)
            ref = self.__get_sol_for_comparison(ref_mem, ud, attribute)

            l2_error = np.linalg.norm(test - ref)
            ref_norm = np.linalg.norm(ref)

            if ref_norm > 0.0:
                rel_l2_error = l2_error / ref_norm
            else:
                rel_l2_error = np.inf if l2_error > self.tolerances[attribute] else 0.0
            max_abs_error = np.max(np.abs(test - ref))

            try:
                assert max_abs_error < self.tolerances[attribute], (
                    "Relative L2 error for attribute %s of %s exceeds tolerance:\n"
                    "L2 error: %.6e\nRelative L2 error: %.6e\nMax abs error: %.6e\nTolerance: %.6e"
                    % (
                        attribute,
                        self.current_run,
                        l2_error,
                        rel_l2_error,
                        max_abs_error,
                        self.tolerances[attribute],
                    )
                )
                logging.info(
                    f"Test passed for {attribute} | "
                    f"L2: {l2_error:.2e}, Rel L2: {rel_l2_error:.2e}, Max Abs: {max_abs_error:.2e}"
                )
            except AssertionError as e:
                logging.info(str(e))
                raise

        logging.info(
            f"""
            {'#' * 10}
            Test passed for {self.current_run}
            {'#' * 10}
            """.strip()
                )

    def __init(self, ds: DiagnosticState):
        tp = test_params(ds)

        self.tps = {
            ds.test_name: tp,
        }
        # self.tps = [tv_2D]

    def __get_tc(self):
        self.tcs = {}
        for test_name, test_param in self.tps.items():
            fn = test_param.fn + ".h5"
            tc = vis_utils.test_case(
                fn, test_param.dir, test_param.Nx, test_param.Ny, ""
            )

            self.tcs[test_name] = tc

    def __read_yaml(self):
        with open("./src/tests/test_targets.yml", "r") as infile:
            self.target = yaml.safe_load(infile)

    def __plot_comparison(self, mem, ref_mem, ud):
        tp = self.tps[self.current_run]

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
            test_sol = np.copy(getattr(Sol, attribute).T)
            if attribute != "rho":
                rho = getattr(Sol, "rho").T
                test_sol /= rho

                # if attribute == 'rhoY':
                #     test_sol -= mpv.HydroState.Y0[:,np.newaxis]

        else:
            # test_sol = mpv.p2_nodes.T * ud.Msq
            # test_sol -= mpv.HydroState_n.pi0[:,np.newaxis]
            test_sol = get_p_from_pressure_related_fields(mem, ud).T
            # pass

        return test_sol

    @staticmethod
    def __get_ens(tc, params, attribute, time_increment=False, summed=True, normed=False):
        if time_increment and attribute == "p2_nodes":
            times = [params.times[0]-1, params.times[0]]
        else:
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
            fn=tc.base_fn,
            load_ic=False,
            avg=False,
        )

        if time_increment and attribute == "p2_nodes":
            ens = ens[1] - ens[0]
        else:
            ens = ens[0] # removes time axis
        ens = ens[0] # removes ensemble axis

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
