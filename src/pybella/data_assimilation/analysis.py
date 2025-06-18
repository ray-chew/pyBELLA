import logging
import numpy as np

from ..utils import sim_params as params
from ..flow_solver.utils.boundary import cell_boundary as bdry_c
from ..flow_solver.utils.boundary import node_boundary as bdry_n

from . import (
    etpf as da_etpf,
    post_processing as da_post_processing,
    letkf as da_letkf,
    utils as da_utils,
)


def do_for_window(tout, outer_step, results, sst, writer):
    # mp = sst.model_params
    dp = sst.da_params
    # ens = dp.sol_ens

    ######################################################
    # Analysis step
    ######################################################
    tout = np.around(tout, 3)
    if sst.N > 1 and tout in dp.dap.da_times:
        futures = []

        ######################################################
        # Update ensemble with forecast
        ######################################################
        for mem in results:
            elem, node, sol, _, npf, th, _ = mem
            bdry_c.set_ghost_cells(mem, sst.ud)
            bdry_n.set_ghost_nodes(npf.p2_nodes, node, sst.ud)

        # ens.set_members(results, tout)
        sst.ensemble_state.set_members(results)

        ######################################################
        # Write output before assimilating data
        ######################################################
        logging.info("Starting output...")
        for mem in sst.ensemble_state:
            elem, node, sol, _, npf, th, _ = mem
            if params.label_type == "STEP":
                step = outer_step
                label = "ensemble_mem=%i_%.3d" % (n, step)
            else:
                label = "ensemble_mem=%i_%.3f" % (n, tout)
            writer.write_all(sol, npf, elem, node, th, str(label) + "_before_da")

        ##################################################
        # LETKF with batch observations
        ##################################################
        if dp.dap.da_type == "batch_obs":
            logging.info("Starting analysis... for batch observations")
            for attr in dp.dap.obs_attributes:
                logging.info("Assimilating %s..." % attr)
                logging.info("Assimilating %s..." % attr)
                # future = client.submit(da_interface, *[s_res,obs_current,dap.inflation_factor,attr,N,ud,dap.loc[attr]])
                future = da_letkf.da_interface(
                    results, dp.dap, dp.obs, attr, tout, sst.N, sst.ud
                )
                futures.append(future)

            # analysis = client.gather(futures)
            analysis = futures
            # analysis = np.array(analysis)

            logging.info("Writing analysis...")
            cnt = 0
            for attr in dp.dap.obs_attributes:
                current = analysis[cnt]
                for n in range(sst.N):
                    setattr(results[:, dp.dap.loc[attr], ...][n], attr, current[n])
                cnt += 1

        ##################################################
        # LETKF with grid-point localisation
        ##################################################
        elif dp.dap.da_type == "rloc":
            logging.info("Starting analysis... for rloc algorithm")
            elem, node = sst.ensemble_state.get_grid()
            results = da_utils.HSprojector_3t2D(results, elem, dp.dap, sst.N)
            results = dp.rloc.analyse(
                results, dp.obs, dp.obs_covar, dp.obs_mask, sst.N, tout
            )
            results = da_utils.HSprojector_2t3D(results, elem, node, dp.dap, sst.N)
            # if hasattr(dap, 'converter'):
            # results = dap.converter(results, N, npf, elem, node, th, ud)

        ##################################################
        # ETPF
        ##################################################
        elif dp.dap.da_type == "etpf":
            da_utils.ensemble_inflation(
                results, dp.dap.attributes, dp.dap.inflation_factor, sst.N
            )
            results = da_etpf.da_interface(
                results,
                dp.obs,
                dp.dap.obs_attributes,
                dp.dap.rejuvenation_factor,
                dp.dap.da_times,
                tout,
                sst.N,
            )

        ##################################################
        # Post-processing
        ##################################################
        elif dp.dap.da_type == "pprocess":
            results = da_post_processing.interface()

        else:
            assert 0, "DA type not implemented: use 'rloc', 'batch_obs' or 'etpf'."

    ######################################################
    # Update ensemble with analysis
    ######################################################
    for mem in results:
        elem, node, sol, npf, th, _, _ = mem
        bdry_c.set_ghost_cells(mem, sst.ud)
        p2_nodes = npf.p2_nodes
        bdry_n.set_ghost_nodes(p2_nodes, node, sst.ud)

    sst.ensemble_state.set_members(results)
