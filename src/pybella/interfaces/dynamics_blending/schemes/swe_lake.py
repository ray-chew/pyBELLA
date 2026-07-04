"""SWE <-> lake blending conversions.

The shallow-water cases run through the gas-dynamics solver in the 2D x-y
plane (gamma = 2 equivalence, see ``tests/test_swe_vortex.py``); the depth
field lives on ``sol.rho`` and the free-surface perturbation on
``npf.p2_nodes``. ``do_swe_to_lake_conv`` freezes the free surface (SWE ->
rigid-lid lake reference depth); ``do_lake_to_swe_conv`` releases the lid
after a look-ahead lake step, blending the nodal pressure by
``ud.blending_weight`` — the shallow-water analogue of the comp <-> psinc
pair in ``comp_psinc.py``.

The legacy thin-y 3D form of these routines (x, y = one-cell vertical
shell, z) predates the ModelState refactor and is preserved in git history
(last in "refactor pt 7"); this module is its 2D x-y port.
"""

import copy
import logging

import numpy as np
from scipy import signal

from ....utils import io


def _node_to_cell(field_n):
    """Average a nodal field onto cells (2x2 'valid' convolution)."""
    kernel = np.ones((2, 2))
    kernel /= kernel.sum()
    return signal.convolve(field_n, kernel, mode="valid")


def do_swe_to_lake_conv(mem, ud, writer, label):
    """Freeze the free surface: SWE depth -> lake (rigid-lid) reference depth.

    Stores the reference depth on ``ud.mean_val``; it must persist across the
    lake window — ``do_lake_to_swe_conv`` consumes it when the lid is
    released.
    """
    logging.info("swe to lake conversion...")
    sol = mem.sol

    H10 = np.copy(mem.npf.p2_nodes)
    H10 -= H10.mean()
    H10 = _node_to_cell(H10)

    setattr(ud, "mean_val", sol.rho - ud.Msq * H10)

    sol.rhou[...] = sol.rhou / sol.rho * ud.mean_val
    sol.rhov[...] = sol.rhov / sol.rho * ud.mean_val
    sol.rhow[...] = sol.rhow / sol.rho * ud.mean_val
    sol.rhoY[...] = sol.rhoY / sol.rho * ud.mean_val
    sol.rho[...] = ud.mean_val


def do_lake_to_swe_conv(mem, ud, label, writer, step, tout):
    """Release the rigid lid: lake -> SWE.

    Runs a look-ahead lake step to ``tout`` (structurally parallel to
    ``do_psinc_to_comp_conv``), blends the nodal pressure between the
    look-ahead result and the frozen pre-step state per
    ``ud.blending_weight`` / ``ud.blending_type``, then reconstructs the SWE
    depth and momenta from the blended pressure around ``ud.mean_val``.

    Like ``do_psinc_to_comp_conv``, the look-ahead clock advance is rolled
    back: the legacy scheduling passed t/step by value, so the outer clock
    stayed untouched.
    """
    from ....flow_solver.discretisation import time_update

    logging.info("doing lake-to-swe time-update...")
    sol_freeze = copy.deepcopy(mem.sol)
    npf_freeze = copy.deepcopy(mem.npf)
    time_freeze = (mem.time.t, mem.time.step, mem.time.window_step)
    # reference clock for the look-ahead ([0, step] in the paper-era
    # data.time_update): window_step = 0 keeps the eos schedule in the
    # limit (lake) regime under continuous blending
    mem.time.window_step = 0

    # exactly ONE look-ahead step — same ULP knife-edge as
    # do_psinc_to_comp_conv (see the note there)
    stepmax_freeze = ud.stepmax
    ud.stepmax = mem.time.step + 1

    ret = time_update.do(
        mem,
        ud,
        tout,
        bld=None,
        writer=None,
        debug_writer=io.NullDebugWriter(),
    )
    ud.stepmax = stepmax_freeze
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

    if writer is not None:
        writer.populate(str(label) + "_after_full_step", "dp2n", dp2n)

    mem.sol = sol_freeze
    mem.npf = npf_freeze
    mem.npf.p2_nodes[...] = dp2n

    logging.info("lake to swe conversion...")
    H10 = np.copy(mem.npf.p2_nodes)
    H10 -= H10.mean()
    H1 = ud.mean_val + ud.Msq * _node_to_cell(H10)

    sol = mem.sol
    sol.rho[...] = H1
    sol.rhou[...] = sol.rhou / ud.mean_val * sol.rho
    sol.rhov[...] = sol.rhov / ud.mean_val * sol.rho
    sol.rhow[...] = sol.rhow / ud.mean_val * sol.rho
    sol.rhoY[...] = sol.rhoY / ud.mean_val * sol.rho

    return mem
