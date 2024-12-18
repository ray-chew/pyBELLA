# some diagnostics
import time
import logging

import numpy as np
import termcolor

# dependencies of the atmospheric flow solver
from .dycore.discretisation import time_update    as dis_time_update
from .dycore.utils import boundary as bdry

# dependencies of the interface subpackage
from .interfaces.dynamics_blending import prepare as blending_prepare

# dependencies of the data assimilation subpackage
from .data_assimilation import (
    prepare as da_prepare,
    etpf as da_etpf,
    post_processing as da_post_processing,
    letkf as da_letkf,
    utils as da_utils
)

# input file
from .utils.sim_params import debug
from .utils import prepare
from .utils import io as io

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
            if N > 1:
                mem[3][0] = 0 if tout_old in dap.da_times else mem[3][0]
            if N == 1:
                mem[3][0] = mem[3][1]
            logging.info("For ensemble member = %i..." % mem_cnt)
            future = dis_time_update.do(
                mem[0],
                mem[1],
                mem[2],
                t,
                tout,
                ud,
                elem,
                node,
                mem[3],
                th,
                blend,
                step_writer,
                debug,
            )

            if ud.diag:
                diag_comparison.test_do(
                    future[0], future[2].p2_nodes, plot=ud.diag_plot_compare
                )

            futures.append(future)
            mem_cnt += 1

        # Dask commands, used only when parallelisation is
        # enabled
        # results = client.gather(futures)
        results = np.copy(futures)
        results = np.array(results)
        # s_res = client.scatter(results)

        ######################################################
        # Analysis step
        ######################################################
        tout = np.around(tout, 3)
        if N > 1 and tout in dap.da_times:
            futures = []

            ######################################################
            # Update ensemble with forecast
            ######################################################
            for n in range(N):
                Sol = results[n][dap.loc_c]
                bdry.set_explicit_boundary_data(Sol, elem, ud, th, mpv)
                results[n][dap.loc_c] = Sol
                p2_nodes = getattr(results[n][dap.loc_n], "p2_nodes")
                bdry.set_ghostnodes_p2(p2_nodes, node, ud)
                setattr(results[n][dap.loc_n], "p2_nodes", p2_nodes)

            ens.set_members(results, tout)

            ######################################################
            # Write output before assimilating data
            ######################################################
            logging.info(termcolor.colored("Starting output...", "yellow"))
            for n in range(N):
                Sol = ens.members(ens)[n][0]
                mpv = ens.members(ens)[n][2]

                if label_type == "STEP":
                    step = outer_step
                    label = "ensemble_mem=%i_%.3d" % (n, step)
                else:
                    label = "ensemble_mem=%i_%.3f" % (n, tout)
                writer.write_all(Sol, mpv, elem, node, th, str(label) + "_before_da")

            ##################################################
            # LETKF with batch observations
            ##################################################
            if dap.da_type == "batch_obs":
                logging.info("Starting analysis... for batch observations")
                for attr in dap.obs_attributes:
                    logging.info("Assimilating %s..." % attr)
                    logging.info("Assimilating %s..." % attr)
                    # future = client.submit(da_interface, *[s_res,obs_current,dap.inflation_factor,attr,N,ud,dap.loc[attr]])
                    future = da_letkf.da_interface(results, dap, obs, attr, tout, N, ud)
                    futures.append(future)

                # analysis = client.gather(futures)
                analysis = futures
                # analysis = np.array(analysis)

                logging.info("Writing analysis...")
                cnt = 0
                for attr in dap.obs_attributes:
                    current = analysis[cnt]
                    for n in range(N):
                        setattr(results[:, dap.loc[attr], ...][n], attr, current[n])
                    cnt += 1

            ##################################################
            # LETKF with grid-point localisation
            ##################################################
            elif dap.da_type == "rloc":
                logging.info(
                    termcolor.colored("Starting analysis... for rloc algorithm", "green")
                )
                results = da_utils.HSprojector_3t2D(results, elem, dap, N)
                results = rloc.analyse(results, obs, obs_covar, obs_mask, N, tout)
                results = da_utils.HSprojector_2t3D(results, elem, node, dap, N)
                # if hasattr(dap, 'converter'):
                # results = dap.converter(results, N, mpv, elem, node, th, ud)

            ##################################################
            # ETPF
            ##################################################
            elif dap.da_type == "etpf":
                da_utils.ensemble_inflation(results, dap.attributes, dap.inflation_factor, N)
                results = da_etpf.da_interface(
                    results,
                    obs,
                    dap.obs_attributes,
                    dap.rejuvenation_factor,
                    dap.da_times,
                    tout,
                    N,
                )

            ##################################################
            # Post-processing
            ##################################################
            elif dap.da_type == "pprocess":
                results = da_post_processing.interface()

            else:
                assert 0, "DA type not implemented: use 'rloc', 'batch_obs' or 'etpf'."

        ######################################################
        # Update ensemble with analysis
        ######################################################
        for n in range(N):
            Sol = results[n][dap.loc_c]
            bdry.set_explicit_boundary_data(Sol, elem, ud, th, mpv)
            results[n][dap.loc_c] = Sol
            p2_nodes = getattr(results[n][dap.loc_n], "p2_nodes")
            bdry.set_ghostnodes_p2(p2_nodes, node, ud)
            setattr(results[n][dap.loc_n], "p2_nodes", p2_nodes)

        ens.set_members(results, tout)

        ######################################################
        # Write output at tout
        ######################################################
        logging.info(termcolor.colored("Starting output...", "yellow"))
        for n in range(N):
            Sol = ens.members(ens)[n][0]
            mpv = ens.members(ens)[n][2]

            if label_type == "STEP":
                step = outer_step
                label = "ensemble_mem=%i_%.3d" % (n, step)
            else:
                label = "ensemble_mem=%i_%.3f" % (n, tout)
            writer.write_all(Sol, mpv, elem, node, th, str(label) + "_after_full_step")

        # synchronise_variables(mpv, Sol, elem, node, ud, th)
        t = tout
        tout_old = np.copy(tout)
        logging.info(termcolor.colored("tout = %.3f" % tout, "yellow"))

        tout_cnt += 1
        outer_step += 1
        if outer_step > ud.stepmax:
            break

    toc = time.time()
    logging.info(termcolor.colored("Time taken = %.6f" % (toc - tic), "yellow"))

    writer.close_everything()


if __name__ == "__main__":
    main()