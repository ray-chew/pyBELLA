import logging
import copy

import numpy as np
from scipy import signal

from ...flow_solver.physics import eos as gd_eos
from ...utils import io


class Blend(object):
    """
    Class that takes care of the blending interface.
    """

    def __init__(self, ud):
        self.bb = False
        self.cb = ud.continuous_blending
        self.psinc_init = ud.no_of_pi_initial
        # self.psinc_trans = ud.no_of_pi_transition
        self.hydro_init = ud.no_of_hy_initial

        if self.psinc_init > 0 and self.cb:
            self.bb = True

        self.c_init = self.criterion_init
        # self.c_trans = self.criterion_trans

        self.fac = ud.Msq

    def criterion_init(self, step):
        return step == (self.psinc_init) and self.cb and self.bb

    # def criterion_trans(self, step):
    #     return step <= self.psinc_trans and self.cb and not self.bb

    def convert_p2n(self, p2n):
        ndim = p2n.ndim
        dp2n = p2n - p2n.mean()

        self.kernel = np.ones([2] * ndim)
        dp2c = signal.fftconvolve(dp2n, self.kernel, mode="valid") / self.kernel.sum()

        # self.dp2n = dp2n - dp2n.mean()
        # self.dp2c = dp2c - dp2c.mean()

        self.dp2n = dp2n
        self.dp2c = dp2c

        return dp2c

    def update_sol(self, mem, ud, sgn, label=None, writer=None):
        if writer != None:
            writer.populate(str(label) + "_before_blending", "dp2n", self.dp2n)
        if writer != None:
            writer.write_all(mem, str(label) + "_before_blending")

        sol = mem.sol
        npf = mem.npf
        th = mem.th

        if sgn == "bef":
            sign = -1.0
        elif sgn == "aft":
            sign = +1.0
        else:
            assert 0, "sgn == bef or sgn == aft"

        rho = np.copy(sol.rho)
        rhoY = np.copy(sol.rhoY)

        Y = rhoY / rho

        if ud.blending_mean == "rhoY":
            rhoYc = (rhoY**th.gm1 + sign * self.fac * self.dp2c) ** (th.gm1inv)
        elif ud.blending_mean == "1.0":
            rhoYc = (1.0 + sign * self.fac * self.dp2c) ** (th.gm1inv)

        alpha = rhoYc / sol.rhoY

        if ud.blending_conv == "rho":
            ### keep theta, convert rho
            sol.rho[...] = rho * alpha
            sol.rhoY[...] = sol.rho * Y

            rho_fac = sol.rho / rho
            sol.rhou[...] *= rho_fac
            sol.rhov[...] *= rho_fac
            sol.rhow[...] *= rho_fac
            sol.rhoX[...] *= rho_fac

        elif ud.blending_conv == "theta":
            ### keep rho, convert theta
            Yc = Y * alpha
            sol.rhoY[...] = rho * Yc
            sol.rhoX[...] = rho * (1.0 / Yc - npf.HydroState.S0.reshape(1, -1))
        else:
            assert 0, "ud.blending_conv undefined."

        if writer != None:
            writer.write_all(mem, str(label) + "_after_blending")

    def update_p2n(self, npf):
        npf.p2_nodes = self.dp2n


######################################################
# COMP - PSINC blending
######################################################


def do_comp_to_psinc_conv(mem, bld, ud, label, writer):
    logging.info("Converting COMP to PSINC")
    dp2n = mem.npf.p2_nodes
    bld.convert_p2n(dp2n)
    bld.update_sol(mem, ud, "bef", label=label, writer=writer)
    bld.update_p2n(mem.npf)

    return mem


def do_psinc_to_comp_conv(
    mem,
    ud,
    bld,
    label,
    writer,
    step,
    tout,
):
    from ...flow_solver.discretisation import time_update

    logging.info(f"Blending... step = {step}")
    sol_freeze = copy.deepcopy(mem.sol)
    npf_freeze = copy.deepcopy(mem.npf)

    ret = time_update.do(
        mem,
        ud,
        tout,
        bld=None,
        writer=None,
        debug_writer=io.NullDebugWriter(),
    )

    fac_old = ud.blending_weight
    fac_new = 1.0 - fac_old
    dp2n_0 = fac_new * ret.npf.p2_nodes_half + fac_old * npf_freeze.p2_nodes_half
    dp2n_1 = fac_new * ret.npf.p2_nodes + fac_old * npf_freeze.p2_nodes

    if ud.blending_type == "half":
        dp2n = dp2n_0
    elif ud.blending_type == "full":
        dp2n = dp2n_1
    else:
        assert 0, "incorrect ud.blending_type"

    if writer != None:
        writer.populate(
            str(label) + "_after_full_step", "p2_start", npf_freeze.p2_nodes
        )
    if writer != None:
        writer.populate(str(label) + "_after_full_step", "p2_end", ret.npf.p2_nodes)
    mem.sol = sol_freeze
    mem.npf = npf_freeze

    # elem, node, _, _, _, th, _ = mem

    if writer != None:
        writer.populate(str(label) + "_after_full_step", "dp2n", dp2n)
    logging.info("Converting PSINC to COMP")
    bld.convert_p2n(dp2n)
    bld.update_sol(mem, ud, "aft", label=label, writer=writer)
    bld.update_p2n(mem.npf)

    # mem.time.step -= 1
    # mem.time.window_step -= 1

    return mem


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
    from ...flow_solver.discretisation import time_update

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


######################################################
# Nonhydrostatic - Hydrostatic blending
######################################################
def do_nonhydro_to_hydro_conv(
    sol, flux, npf, bld, elem, node, th, ud, label, writer, step, window_step, t, dt
):
    logging.info("nonhydrostatic to hydrostatic conversion...")
    # bld.convert_p2n(npf.p2_nodes)
    # bld.update_sol(sol,elem,node,th,ud,npf,'bef',label=label,writer=writer)
    # sol.rhov = sol.rhov_half

    # sol_tmp = deepcopy(sol)
    # flux_tmp = deepcopy(flux)
    # npf_tmp = deepcopy(npf)

    # nonhydro to hydro blending incomplete.
    # ret = data.time_update(sol,flux,npf, t, t+1*dt, ud, elem, node, [0,0], th, bld=None, writer=None, debug=False)

    # sol = sol_tmp
    # flux = flux_tmp
    # npf = npf_tmp
    # sol = ret[0]
    # flux = ret[1]
    # npf = ret[2]
    # sol = deepcopy(ret[0])
    # npf = deepcopy(ret[2])
    # sol.rhov[...] = sol.rhov_half
    # t += 0.5*dt
    # t += 1*dt
    return sol, npf, t


def do_hydro_to_nonhydro_conv(
    sol, flux, npf, bld, elem, node, th, ud, label, writer, step, window_step, t, dt
):
    logging.info("hydrostatic to nonhydrostatic conversion...")
    logging.info(f"Blending... step = {step}")

    # sol_tmp = deepcopy(sol)
    # flux_tmp = deepcopy(flux)
    # npf_tmp = deepcopy(npf)

    # ret = data.time_update(sol,flux,npf, t, t+dt, ud, elem, node, [0,step-1], th, bld=None, writer=None, debug=False)

    # sol = sol_tmp
    # flux = flux_tmp
    # npf = npf_tmp

    # retv_half = ret[0].rhov_half / ret[0].rho_half
    # retv_full = ret[0].rhov / ret[0].rho

    # solv_half = sol.rhov_half / sol.rho_half
    # solv_full = sol.rhov / sol.rho

    # fac_full = 0.5
    # fac_half = 1.0 - fac_full

    # # logging.info(np.sum(solv_full))
    # # logging.info(np.sum(retv_half))
    # # logging.info(np.sum((fac_full * solv_full + fac_half * retv_half)))
    # # logging.info(np.sum(fac_half * retv_half))

    # fac_full = 0.5
    # fac_half = 0.5

    # # logging.info(np.sum(solv_full))
    # # logging.info(np.sum(retv_half))
    # # logging.info(np.sum((fac_full * solv_full + fac_half * retv_half)))
    # # logging.info(np.sum(fac_half * retv_half))

    # if writer != None: writer.populate(str(label)+'_after_full_step', 'ret_half', ret[0].rhov_half)
    # if writer != None: writer.populate(str(label)+'_after_full_step', 'ret_full', ret[0].rhov)

    # if writer != None: writer.populate(str(label)+'_after_full_step', 'solv_half', sol.rhov_half)
    # if writer != None: writer.populate(str(label)+'_after_full_step', 'solv_full', sol.rhov)

    # sol.rhov = sol.rho * (fac_full * solv_full + fac_half * retv_half)
    # if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_end', ret[2].p2_nodes)

    # fac_npf_half = 0.5
    # fac_npf_full = 1.0 - fac_npf_half
    # npf.p2_nodes = fac_npf_half * npf.p2_nodes + fac_npf_full * ret[2].p2_nodes
    # dp2n = ret[2].p2_nodes_half
    # bld.convert_p2n(dp2n)
    # bld.update_sol(sol,elem,node,th,ud,npf,'aft',label=label,writer=writer)
    # bld.update_p2n(sol,npf,node,th,ud)
    #

    ###############################
    # alternative version
    ###############################

    # if c1 or c2:
    #     logging.info(
    #         termcolor.colored("hydrostatic to nonhydrostatic conversion...", "blue")
    #     )

    # writer.write_all(mem, str(label) + "_half_full")
    # writer.populate(str(label) + "_ic", "pwchi", sol.pwchi)

    # if test_hydrob == False:
    #     sol = copy.deepcopy(sol_half_old)
    #     # npf = copy.deepcopy(npf_half_old)

    #     logging.info(termcolor.colored("test_hydrob == False", "red"))
    #     writer.write_all(mem, str(label) + "_quarter")

    #     writer.populate(str(label) + "_quarter", "pwchi", sol.pwchi)

    #     logging.info("quarter dt = %.8f" % (dt * 0.5))

    #     ret = do(
    #         sol_half_old,
    #         flux_half_old,
    #         npf_half_old,
    #         dt - 0.5 * dt,
    #         dt + 0.5 * dt,
    #         ud,
    #         elem,
    #         node,
    #         [0, 0],
    #         th,
    #         bld=None,
    #         writer=None,
    #         debug=False,
    #     )

    #     sol_tu = copy.deepcopy(ret[0])
    #     # npf_tu = copy.deepcopy(ret[2])
    #     sol.rho[...] = sol_tu.rho_half
    #     sol.rhou[...] = sol_tu.rhou_half
    #     sol.rhov[...] = sol_tu.rhov_half
    #     sol.rhow[...] = sol_tu.rhow_half
    #     sol.rhoX[...] = sol_tu.rhoX_half
    #     sol.rhoY[...] = sol_tu.rhoY_half
    #     sol.pwchi[...] = sol_tu.pwchi

    #     # npf.p2_nodes[...] = npf_tu.p2_nodes_half

    #     writer.write_all(mem, str(label) + "_half")

    #     writer.populate(str(label) + "_half", "pwchi", sol.pwchi)

    #     ret = do(
    #         sol,
    #         flux,
    #         npf,
    #         dt,
    #         2.0 * dt,
    #         ud,
    #         elem,
    #         node,
    #         [0, 0],
    #         th,
    #         bld=None,
    #         writer=None,
    #         debug=False,
    #     )

    #     sol = copy.deepcopy(ret[0])
    #     flux = copy.deepcopy(ret[1])
    #     npf = copy.deepcopy(ret[2])

    # if test_hydrob == True:
    #     sol = copy.deepcopy(sol_half_old)
    #     # npf = copy.deepcopy(npf_half_old)

    #     logging.info(termcolor.colored("test_hydrob == False", "red"))
    #     writer.write_all(mem, str(label) + "_quarter")

    #     # writer.populate(str(label)+'_quarter', 'pwchi', sol.pwchi)

    #     logging.info("quarter dt = %.8f" % (dt * 0.5))

    #     ret = do(
    #         sol_half_old,
    #         flux_half_old,
    #         npf_half_old,
    #         dt - 0.5 * dt,
    #         dt + 0.5 * dt,
    #         ud,
    #         elem,
    #         node,
    #         [0, 0],
    #         th,
    #         bld=None,
    #         writer=None,
    #         debug=False,
    #     )

    #     sol_tu = copy.deepcopy(ret[0])
    #     # npf_tu = copy.deepcopy(ret[2])
    #     sol.rho[...] = sol_tu.rho_half
    #     sol.rhou[...] = sol_tu.rhou_half
    #     sol.rhov[...] = sol_tu.rhov_half
    #     sol.rhow[...] = sol_tu.rhow_half
    #     sol.rhoX[...] = sol_tu.rhoX_half
    #     sol.rhoY[...] = sol_tu.rhoY_half
    #     sol.pwchi[...] = sol_tu.pwchi

    #     # npf.p2_nodes[...] = npf_tu.p2_nodes_half

    #     # writer.write_all(sol,npf,elem,node,th,str(label)+'_half')

    #     # writer.populate(str(label)+'_half', 'pwchi', sol.pwchi)

    #     ret = do(
    #         sol,
    #         flux,
    #         npf,
    #         dt,
    #         2.0 * dt,
    #         ud,
    #         elem,
    #         node,
    #         [0, 0],
    #         th,
    #         bld=None,
    #         writer=None,
    #         debug=False,
    #     )

    #     sol = copy.deepcopy(ret[0])
    #     flux = copy.deepcopy(ret[1])
    #     npf = copy.deepcopy(ret[2])
    #     # writer.write_all(sol,npf,elem,node,th,str(label)+'_half')
    #     # writer.populate(str(label)+'_half', 'pwchi', sol.pwchi)

    #     logging.info(termcolor.colored("test_hydrob == True", "red"))

    # if test_hydrob == False:
    #     dt *= 2.0
    # if c2:
    # ud.is_nonhydrostatic = 1

    return sol, npf


######################################################
# Blending calls from data.py
######################################################
def blending_before_timestep(
    mem,
    ud,
    bld,
    label,
    writer,
    step,
    window_step,
    t,
    dt,
    swe_to_lake,
    debug,
):
    ######################################################
    # Blending : Do full regime to limit regime conversion
    ######################################################
    # do unpacking
    elem, node, sol, npf, th, _, _ = mem

    # these make sure that we are the correct window step
    if bld is not None and window_step == 0:
        # these make sure that blending switches are on
        if (bld.bb or bld.cb) and ud.blending_conv is not None:
            # these distinguish between SWE and Euler blending
            if ud.blending_conv == "swe":
                do_swe_to_lake_conv(sol, npf, elem, node, ud, th, writer, label, debug)
                swe_to_lake = True
            else:
                mem = do_comp_to_psinc_conv(mem, bld, ud, label, writer)

    ######################################################
    # Blending : Do full steps or transition steps?
    ######################################################
    if bld is not None:
        c_init = bld.criterion_init(window_step)
    else:
        c_init = False

    ######################################################
    # Blending : If full blending steps...
    ######################################################
    # check that blending switches are on
    if c_init and bld.cb and ud.blending_conv is not None:
        # distinguish between Euler and SWE blending
        if ud.blending_conv != "swe":
            do_psinc_to_comp_conv(
                mem,
                ud,
                bld,
                ud,
                label,
                writer,
                step,
                t + dt,
            )

    ######################################################
    # Initial Blending
    ######################################################
    # Is initial blending switch on, and if yes, are we in the 0th time-step?
    if ud.initial_blending == True and step < 1 and bld is not None:
        # Distinguish between SWE and Euler blendings
        if ud.blending_conv != "swe":
            if bld.psinc_init > 0:
                ud.is_compressible = 0
                ud.compressibility = 0.0
                mem = do_comp_to_psinc_conv(mem, bld, ud, label, writer)
            elif bld.hydro_init > 0:
                sol, npf, t = do_nonhydro_to_hydro_conv(
                    sol,
                    flux,
                    npf,
                    bld,
                    elem,
                    node,
                    th,
                    ud,
                    label,
                    writer,
                    step,
                    window_step,
                    t,
                    dt,
                )
                ud.is_nonhydrostatic = 0
                ud.nonhydrostasy = 0.0
        else:
            do_swe_to_lake_conv(sol, npf, elem, node, ud, th, writer, label, debug)
            swe_to_lake = True
            ud.is_compressible = 0
            ud.compressibility = 0.0

    # Elif, is initial blending switch on and are we on the 1st time-step?
    # If we are on the first time-step, do we do comp-psinc blending?
    elif (
        ud.initial_blending == True and step == ud.no_of_pi_initial and bld is not None
    ):
        # Distinguish between SWE and Euler blendings
        if ud.blending_conv != "swe":
            do_psinc_to_comp_conv(
                mem,
                ud,
                bld,
                label,
                writer,
                step,
                t + dt,
            )
            ud.is_compressible = 1
            ud.compressibility = 1.0
    # Else, do we do nonhydrostatic-hydrostatic blending?
    elif (
        ud.initial_blending == True and step == ud.no_of_hy_initial and bld is not None
    ):
        if ud.blending_conv != "swe":
            sol, npf = do_hydro_to_nonhydro_conv(
                sol,
                flux,
                npf,
                bld,
                elem,
                node,
                th,
                ud,
                label,
                writer,
                step,
                window_step,
                t,
                dt,
            )
            ud.is_nonhydrostatic = 1
            ud.nonhydrostasy = 1.0
    else:
        ud.is_compressible = gd_eos.is_compressible(ud, window_step)
        ud.compressibility = gd_eos.compressibility(ud, t, window_step)
        ud.is_nonhydrostatic = gd_eos.is_nonhydrostatic(ud, window_step)
        ud.nonhydrostasy = gd_eos.nonhydrostasy(ud, t, window_step)

    return swe_to_lake, sol, npf, t


def blending_after_timestep(
    sol,
    npf,
    bld,
    elem,
    node,
    th,
    ud,
    label,
    writer,
    step,
    window_step,
    t,
    dt,
    swe_to_lake,
    debug,
):
    ######################################################
    # Blending : Are we in the lake regime? And is this
    #            the window step where we go back to SWE?
    ######################################################
    if bld is not None and swe_to_lake and step > 0:
        initialise_lake_to_swe_conv = bld.criterion_init(window_step + 1)

    ######################################################
    # Blending : If we are in the lake regime, is blending
    #            on? If yes, do lake-to-swe conversion.
    ######################################################
    if (
        ud.blending_conv == "swe"
        and swe_to_lake
        and initialise_lake_to_swe_conv
        and bld is not None
    ):
        tmp_CFL = np.copy(ud.CFL)
        ud.CFL = 0.8
        sol, npf = do_lake_to_swe_conv(
            sol,
            flux,
            npf,
            elem,
            node,
            ud,
            th,
            writer,
            label,
            debug,
            step,
            window_step,
            t,
            dt,
        )
        ud.CFL = tmp_CFL
        ud.is_compressible = 1
        ud.compressibility = 1.0

    return sol, npf


def prepare_blending(
    mem,
    ud,
    bld,
    label,
    writer,
    step,
    window_step,
    t,
    dt,
    swe_to_lake,
    debug,
):

    if check_and_apply_initial_hydrostatic_conversion(step, ud, bld):
        ud.is_nonhydrostatic = 0

    swe_to_lake, sol, npf, t = blending_before_timestep(
        mem,
        ud,
        bld,
        label,
        writer,
        step,
        window_step,
        t,
        dt,
        swe_to_lake,
        debug,
    )

    return swe_to_lake, sol, npf, t


def check_and_apply_initial_hydrostatic_conversion(step, ud, bld):
    """
    Check and apply initial blending conversion if needed.

    Returns:
        bool: True if conversion was applied, False otherwise
    """
    if step != 0 or bld is None or "imbal" not in ud.aux:
        return False

    hydrostatic_case = ud.is_nonhydrostatic == 0
    nonhydrostatic_case = ud.is_nonhydrostatic == 1 and ud.initial_blending == True

    if hydrostatic_case or nonhydrostatic_case:
        logging.info("nonhydrostatic to hydrostatic conversion...")
        ud.is_nonhydrostatic = 0
        return True

    return False
