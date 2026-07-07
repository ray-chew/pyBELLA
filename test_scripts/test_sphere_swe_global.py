"""Stage F, increment F5: Williamson TC2 on the FULL pole-to-pole sphere.

The acid test that the pole ghost exchange (F1), pole-face flux closure
(F2), polar filter (F3) and elliptic pole collapse (F4) compose in a real
forecast. TC2's steady zonal flow vanishes at the poles and is
longitude-independent, so the exact state is pole-regular and the filter is
a no-op — any drift is our discretisation.

Short tripwires (a coarse resolution keeps CI fast; the 12-day l2 <= 1e-3
Williamson validation runs from run_scripts). Gates:

- steady geostrophic balance: l2(h) drift stays at the adjustment/Krylov
  floor over 10 steps (calibrated well below the ~5e-4 a wrong Coriolis
  sign/factor produces, as in the channel case);
- the resting global shell is an exact discrete steady state;
- the tangent-plane constraint holds to machine precision (the filter does
  not break it: surface_constraint is re-applied after filtering).
"""

import logging

import numpy as np

logging.disable(logging.INFO)

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.tests import test_sphere_swe_tc2_global as tc2g
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.io.debug import NullDebugWriter


def _run(nsteps, rest=False, nlam=48, nphi=54):
    udo = tc2g.UserData()
    udo.inx, udo.inz = nlam + 1, nphi + 1  # coarse for the tripwire
    udo.stepmax = nsteps
    udo.diag = False
    udo.output_timesteps = False
    udo.output_suffix = "_%i_%i" % (udo.inx - 1, udo.inz - 1)
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = tc2g.sol_init(sol, npf, elem, node, th, ud)
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
    mem = time_update.do(
        mem, ud, tout=1.0e9, bld=None, writer=None, debug_writer=NullDebugWriter()
    )
    assert mem.time.step == nsteps
    return mem, inner, h0


def test_tc2_global_steady_and_tangency():
    mem, inner, h0 = _run(10)
    assert np.isfinite(mem.sol.rho).all()
    dh = mem.sol.rho[inner] - h0
    l2 = np.sqrt(np.mean(dh**2)) / np.sqrt(np.mean(h0**2))
    # measured ~6e-6 (adjustment/Krylov floor); a wrong Coriolis sign/factor
    # drifts ~5e-4, ~100x the gate
    assert l2 < 5.0e-5, l2

    moms = (mem.sol.rhou, mem.sol.rhov, mem.sol.rhow)
    e = mem.elem.metric.e_up
    tangency = np.max(np.abs(sum(moms[k] * e[k] for k in range(3))[inner]))
    assert tangency < 1.0e-11, tangency


def test_tc2_global_rest_state_is_exact():
    mem, inner, h0 = _run(10, rest=True)
    dh = mem.sol.rho[inner] - h0
    assert np.max(np.abs(dh)) < 1.0e-11, np.max(np.abs(dh))
    moms = (mem.sol.rhou, mem.sol.rhov, mem.sol.rhow)
    assert max(np.max(np.abs(m[inner])) for m in moms) < 1.0e-11
