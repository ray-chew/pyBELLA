import logging
import copy

from ....utils import io

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
    from ....flow_solver.discretisation import time_update

    logging.info(f"Blending... step = {step}")
    sol_freeze = copy.deepcopy(mem.sol)
    npf_freeze = copy.deepcopy(mem.npf)
    # the pressure-extraction step below is a throwaway: it must not advance
    # the real clock (pre-ModelState code passed t/step by value; swe_lake.py
    # applies the same freeze/restore)
    time_freeze = (mem.time.t, mem.time.step, mem.time.window_step)
    # run it on the reference's clock ([0, step-1] in the paper-era
    # data.time_update): window_step = 0 keeps the eos schedule in the limit
    # regime — with the live window_step == no_of_pi_initial, continuous
    # blending would flip the throwaway step to compressible and extract a
    # compressible (unprojected) half-time pressure
    mem.time.step -= 1
    mem.time.window_step = 0

    ret = time_update.do(
        mem,
        ud,
        tout,
        bld=None,
        writer=None,
        debug_writer=io.NullDebugWriter(),
    )
    mem.time.t, mem.time.step, mem.time.window_step = time_freeze

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

    if writer != None:
        writer.populate(str(label) + "_after_full_step", "dp2n", dp2n)
    logging.info("Converting PSINC to COMP")
    bld.convert_p2n(dp2n)
    bld.update_sol(mem, ud, "aft", label=label, writer=writer)
    bld.update_p2n(mem.npf)

    return mem
