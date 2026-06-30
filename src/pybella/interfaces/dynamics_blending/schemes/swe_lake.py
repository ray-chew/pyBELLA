import logging
import copy

import numpy as np
from scipy import signal

######################################################
# SWE - Lake blending
######################################################


def do_swe_to_lake_conv(sol, npf, elem, node, ud, th, writer, label, debug):
    logging.info("swe to lake conversion...")

    H1 = sol.rho[
        :,
        2:-2:,
    ][:, 0, :]
    # setattr(ud,'mean_val',H1.mean())

    H10 = npf.p2_nodes[:, 2:-2, :].mean(axis=1)
    H10 -= H10.mean()

    # define 2D kernel
    kernel = np.ones((2, 2))
    kernel /= kernel.sum()

    # do node-to-cell averaging
    H10 = signal.convolve(H10, kernel, mode="valid")

    # H1 = (H1 - ud.mean_val)
    H1 = H1 - ud.Msq * H10
    H1 = np.expand_dims(H1, axis=1)
    H1 = np.repeat(H1, elem.icy, axis=1)
    setattr(ud, "mean_val", H1)

    sol.rhou[...] = sol.rhou / sol.rho * ud.mean_val
    sol.rhov[...] = sol.rhov / sol.rho * ud.mean_val
    sol.rhow[...] = sol.rhow / sol.rho * ud.mean_val
    sol.rhoY[...] = sol.rhoY / sol.rho * ud.mean_val
    sol.rho[...] = ud.mean_val

    # boundary.set_ghostnodes_p2(npf.p2_nodes,node,ud)

    if debug == True:
        writer.write_all(sol, npf, elem, node, th, str(label) + "_after_swe_to_lake")


def do_lake_to_swe_conv(
    sol, flux, npf, elem, node, ud, th, writer, label, debug, step, window_step, t, dt
):
    from ....flow_solver.discretisation import time_update

    if debug == True:
        writer.write_all(sol, npf, elem, node, th, str(label) + "_after_lake_time_step")

    sol_freeze = copy.deepcopy(sol)
    npf_freeze = copy.deepcopy(npf)

    logging.info("doing lake-to-swe time-update...")
    ret = time_update.time_update(
        sol,
        flux,
        npf,
        t,
        t + dt,
        ud,
        elem,
        node,
        [0, step],
        th,
        bld=None,
        writer=None,
        debug=False,
    )

    fac_old = ud.blending_weight
    fac_new = 1.0 - fac_old

    dp2n_0 = fac_new * ret[2].p2_nodes_half + fac_old * npf_freeze.p2_nodes_half
    dp2n_1 = fac_new * ret[2].p2_nodes + fac_old * npf_freeze.p2_nodes

    if ud.blending_type == "half":
        dp2n = dp2n_0
    elif ud.blending_type == "full":
        dp2n = dp2n_1
    else:
        assert 0, "incorrect ud.blending_type"

    sol = copy.deepcopy(sol_freeze)
    npf = copy.deepcopy(npf_freeze)

    npf.p2_nodes[...] = dp2n

    H10 = npf.p2_nodes[:, 2:-2, :].mean(axis=1)
    logging.info("lake to swe conversion...")
    H10 -= H10.mean()

    # define 2D kernel
    kernel = np.ones((2, 2))
    kernel /= kernel.sum()

    # do node-to-cell averaging
    H1 = signal.convolve(H10, kernel, mode="valid")
    # H1 = ud.mean_val + ud.Msq * H1
    # logging.info(colored(H1.max(), 'red'))

    # project H1 back to horizontal slice with ghost cells
    H1 = np.expand_dims(H1, axis=1)
    H1 = np.repeat(H1, elem.icy, axis=1)
    H1 = ud.mean_val + ud.Msq * H1

    sol.rho[...] = H1
    sol.rhou[...] = sol.rhou / ud.mean_val * sol.rho
    sol.rhov[...] = sol.rhov / ud.mean_val * sol.rho
    sol.rhow[...] = sol.rhow / ud.mean_val * sol.rho
    sol.rhoY[...] = sol.rhoY / ud.mean_val * sol.rho

    if debug == True:
        writer.write_all(sol, npf, elem, node, th, str(label) + "_after_lake_to_swe")
    return sol, npf
