import logging

import numpy as np

from . import(
    etpf as da_etpf,
    post_processing as da_post_processing,
    letkf as da_letkf,
    utils as da_utils
)

from ..dycore.utils import boundary as bdry

from ..utils import sim_params as params

def do_for_window(tout, outer_step, results, sst, writer):
    mp = sst.model_params
    dp = sst.da_params
    ens = dp.sol_ens

    ######################################################
    # Analysis step
    ######################################################
    tout = np.around(tout, 3)
    if sst.N > 1 and tout in dp.dap.da_times:
        futures = []

        ######################################################
        # Update ensemble with forecast
        ######################################################
        for n in range(sst.N):
            Sol = results[n][dp.dap.loc_c]
            bdry.set_explicit_boundary_data(Sol, mp.elem, mp.ud, mp.th, mpv)
            results[n][dp.dap.loc_c] = Sol
            p2_nodes = getattr(results[n][dp.dap.loc_n], "p2_nodes")
            bdry.set_ghostnodes_p2(p2_nodes, mp.node, sst.ud)
            setattr(results[n][dp.dap.loc_n], "p2_nodes", p2_nodes)

        ens.set_members(results, tout)

        ######################################################
        # Write output before assimilating data
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
            writer.write_all(Sol, mpv, mp.elem, mp.node, mp.th, str(label) + "_before_da")

        ##################################################
        # LETKF with batch observations
        ##################################################
        if dp.dap.da_type == "batch_obs":
            logging.info("Starting analysis... for batch observations")
            for attr in dp.dap.obs_attributes:
                logging.info("Assimilating %s..." % attr)
                logging.info("Assimilating %s..." % attr)
                # future = client.submit(da_interface, *[s_res,obs_current,dap.inflation_factor,attr,N,ud,dap.loc[attr]])
                future = da_letkf.da_interface(results, dp.dap, dp.obs, attr, tout, sst.N, sst.ud)
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
            logging.info(
                "Starting analysis... for rloc algorithm"
            )
            results = da_utils.HSprojector_3t2D(results, mp.elem, dp.dap, sst.N)
            results = dp.rloc.analyse(results, dp.obs, dp.obs_covar, dp.obs_mask, sst.N, tout)
            results = da_utils.HSprojector_2t3D(results, mp.elem, mp.node, mp.dap, sst.N)
            # if hasattr(dap, 'converter'):
            # results = dap.converter(results, N, mpv, elem, node, th, ud)

        ##################################################
        # ETPF
        ##################################################
        elif dp.dap.da_type == "etpf":
            da_utils.ensemble_inflation(results, dp.dap.attributes, dp.dap.inflation_factor, sst.N)
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
    for n in range(sst.N):
        Sol = results[n][dp.dap.loc_c]
        bdry.set_explicit_boundary_data(Sol, mp.elem, sst.ud, mp.th, mp.mpv)
        results[n][dp.dap.loc_c] = Sol
        p2_nodes = getattr(results[n][dp.dap.loc_n], "p2_nodes")
        bdry.set_ghostnodes_p2(p2_nodes, mp.node, sst.ud)
        setattr(results[n][dp.dap.loc_n], "p2_nodes", p2_nodes)

    ens.set_members(results, tout)