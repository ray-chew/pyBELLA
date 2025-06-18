# some diagnostics
import time
import logging

import numpy as np

from .flow_solver.utils import prepare

# dependencies of the atmospheric flow solver
from .flow_solver.discretisation import time_update as dis_time_update

# dependencies of the interface subpackage
from .interfaces.dynamics_blending import prepare as blending_prepare
from .interfaces.postprocessing import strip_target_file as strip_target


# dependencies of the data assimilation subpackage
from .data_assimilation import prepare as da_prepare, analysis as da_analysis

# package imports
from .utils import io, sim_params as params


##########################################################
# Start main looping
##########################################################
def main():
    sim_state = prepare.initialise()
    blending_prepare.initialise(sim_state)
    da_prepare.initialise(sim_state)

    if sim_state.restart:
        prepare.overwrite_init_with_restart(sim_state)

    writer, step_writer = io.initialise(sim_state)

    tic = time.time()

    ######################################################
    # Time looping over data assimilation windows
    ######################################################
    tout_old = -np.inf
    tout_cnt = 0
    outer_step = 0
    for tout in sim_state.ud.tout:
        sst = sim_state
        es = sst.ensemble_state
        dp = sst.da_params

        futures = []

        # based on the initial blending parameter, define if we want to blend for this assimilation window or simulaiton run.
        blend = blending_prepare.init_da_window(sim_state, tout_old, outer_step)

        ######################################################
        # Forecast step
        ######################################################
        logging.info("##############################################")
        logging.info("Next tout = %.3f" % tout)
        logging.info("Starting forecast...")
        for cnt, mem in enumerate(es):
            # handling of DA window step counter
            if sst.N > 1:
                if tout_old in dp.dap.da_times:
                    mem.time.window_step = 0
            if sst.N == 1:
                mem.time.window_step = mem.time.step

            debug_writer = io.create_debug_writer(params.debug, writer, mem)

            logging.info("For ensemble member = %i..." % cnt)
            mem = dis_time_update.do(
                mem, sst.ud, tout, blend, step_writer, debug_writer
            )

            if sst.ud.diag:
                if sst.ud.diag_updt_targets:
                    strip_target.do(
                        step_writer.OUTPUT_FILENAME
                        + step_writer.BASE_NAME
                        + step_writer.SUFFIX
                        + ".h5",
                        sst.ud,
                        time_increment=sst.diag_comparison.time_increment,
                    )
                    sst.diag_comparison.update_targets()
                else:
                    sst.diag_comparison.test_do(mem, sst.ud)

            futures.append(mem)

        # Dask commands, used only when parallelisation is
        # enabled
        results = np.array(futures)

        da_analysis.do_for_window(tout, outer_step, results, sst, writer)

        ######################################################
        # Write output at tout
        ######################################################
        logging.info("Starting output...")
        for n, mem in enumerate(es):
            if params.label_type == "STEP":
                step = outer_step
                label = "ensemble_mem=%i_%.3d" % (n, step)
            else:
                label = "ensemble_mem=%i_%.3f" % (n, tout)
            writer.write_all(mem, str(label) + "_after_full_step")

        tout_old = np.copy(tout)
        logging.info("tout = %.3f" % tout)

        tout_cnt += 1
        outer_step += 1
        if outer_step > sst.ud.stepmax:
            break

    toc = time.time()
    logging.info("Time taken = %.6f" % (toc - tic))

    writer.close_everything()


if __name__ == "__main__":
    main()
