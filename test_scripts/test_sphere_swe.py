"""Physics gates for shallow water on the spherical channel (sphere pt 3).

Short-horizon TRIPWIRES (the golden masters in ``test_flow_solver.py``
gate regression; the 12-day Williamson TC2 l2 <= 1e-3 validation runs
from ``run_scripts`` where wall-clock allows):

- TC2 discrete geostrophic balance: after 10 steps the depth error
  stays at the adjustment-transient floor (~1e-7), THREE orders below
  the imbalance level a wrong Coriolis sign/factor produces (~1e-4 by
  step 3, 2e-3 by step 20 — measured; this pins the rotation-vector
  convention of ``SphericalShellMap.rotation_axis_cart``);
- tangent-plane constraint: max |m . e_r| stays at machine zero;
- the resting shell is an EXACT discrete steady state.
"""

import logging

import numpy as np

logging.disable(logging.INFO)

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.tests import test_sphere_swe_tc2 as tc2
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.io.debug import NullDebugWriter


def _run(nsteps, rest=False):
    udo = tc2.UserData()
    udo.stepmax = nsteps
    udo.diag = False
    udo.output_timesteps = False
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = tc2.sol_init(sol, npf, elem, node, th, ud)
    if rest:
        sol.rho[...] = 1.0
        sol.rhou[...] = 0.0
        sol.rhov[...] = 0.0
        sol.rhow[...] = 0.0
        sol.rhoY[...] = (ud.g_swe / 2.0) ** th.gamminv
        npf.p2_nodes[...] = (ud.g_swe / 2.0) ** th.gamminv
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())

    inner = (slice(2, -2), slice(2, -2), slice(2, -2))
    h0 = sol.rho[inner].copy()
    m0 = [sol.rhou[inner].copy(), sol.rhov[inner].copy(), sol.rhow[inner].copy()]
    mem = time_update.do(
        mem, ud, tout=1.0e9, bld=None, writer=None, debug_writer=NullDebugWriter()
    )
    assert mem.time.step == nsteps
    return mem, inner, h0, m0


def test_tc2_short_balance_and_tangency():
    mem, inner, h0, _ = _run(10)
    dh = mem.sol.rho[inner] - h0
    l2 = np.sqrt(np.mean(dh**2)) / np.sqrt(np.mean(h0**2))
    assert l2 < 1.0e-6, l2  # measured ~5e-8; wrong Coriolis sign: ~5e-4

    moms = (mem.sol.rhou, mem.sol.rhov, mem.sol.rhow)
    e = mem.elem.metric.e_up
    tangency = np.max(np.abs(sum(moms[k] * e[k] for k in range(3))[inner]))
    assert tangency < 1.0e-12, tangency


def test_rest_state_is_exact():
    mem, inner, h0, m0 = _run(10, rest=True)
    assert np.array_equal(mem.sol.rho[inner], h0)
    moms = (mem.sol.rhou, mem.sol.rhov, mem.sol.rhow)
    for k in range(3):
        assert np.max(np.abs(moms[k][inner] - m0[k])) < 1.0e-13
