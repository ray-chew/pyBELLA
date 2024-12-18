# some diagnostics
import time
import logging

import numpy as np

# dependencies of the atmospheric flow solver
from .dycore.discretisation import time_update    as dis_time_update

# dependencies of the interface subpackage
from .interfaces.dynamics_blending import prepare as blending_prepare

# dependencies of the data assimilation subpackage
from .data_assimilation import (
    prepare as da_prepare,
    analysis as da_analysis
)

# package imports
from .utils import (
    prepare,
    io,
    sim_params as params
)

##########################################################
# Start main looping
##########################################################
def main():
    sim_state = prepare.initialise()

    blending_prepare.initialise(sim_state)
    da_prepare.initialise(sim_state)
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
        mp = sst.model_params
        dp = sst.da_params
        ens = dp.sol_ens

        futures = []

        blend = blending_prepare.init_da_window(sim_state, tout_old, outer_step)

        ######################################################
        # Forecast step
        ######################################################
        logging.info("##############################################")
        logging.info("Next tout = %.3f" % tout)
        logging.info("Starting forecast...")
        mem_cnt = 0
        for mem in ens.members(ens):

            # handling of DA window step counter
            if sst.N > 1:
                mem[3][0] = 0 if tout_old in dp.dap.da_times else mem[3][0]
            if sst.N == 1:
                mem[3][0] = mem[3][1]
            logging.info("For ensemble member = %i..." % mem_cnt)
            future = dis_time_update.do(
                mem[0],
                mem[1],
                mem[2],
                sst.t,
                tout,
                sst.ud,
                mp.elem,
                mp.node,
                mem[3],
                mp.th,
                blend,
                step_writer,
                params.debug,
            )

            if sst.ud.diag:
                sst.diag_comparison.test_do(
                    future[0], future[2].p2_nodes, plot=sst.ud.diag_plot_compare
                )

            futures.append(future)
            mem_cnt += 1

        # Dask commands, used only when parallelisation is
        # enabled
        # results = client.gather(futures)
        results = np.copy(futures)
        results = np.array(results)
        # s_res = client.scatter(results)

        da_analysis.do_for_window(tout, outer_step, results, sst, writer)

        ######################################################
        # Write output at tout
        ######################################################
        logging.info("Starting output...")
        for n in range(sst.N):
            Sol = ens.members(ens)[n][0]
            mpv = ens.members(ens)[n][2]

            if params.label_type == "STEP":
                step = outer_step
                label = "ensemble_mem=%i_%.3d" % (n, step)
            else:
                label = "ensemble_mem=%i_%.3f" % (n, tout)
            writer.write_all(Sol, mpv, mp.elem, mp.node, mp.th, str(label) + "_after_full_step")

        # synchronise_variables(mpv, Sol, elem, node, ud, th)
        t = tout
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