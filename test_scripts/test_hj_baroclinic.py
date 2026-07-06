"""Hughes & Jablonowski (2023) pt 1 — well-balanced-background steadiness gate.

The Ullrich base state (``pybella.tests.test_hj_baroclinic``) mapped onto the
deep spherical shell WITHOUT topography must stay (nearly) steady: gradient-
wind balance is a steady state of the continuous equations, so the only
motion the discrete solver should generate over a short run is the small
truncation-level adjustment of the balance — NOT a growing meridional
circulation and certainly not a blow-up. This is the 3D analogue of the TC2
balance tripwire: a wrong Coriolis sign / factor or a mis-placed meridional
pressure gradient breaks geostrophy and shows up immediately as spurious
meridional wind.

The gate quantifies three things over a coarse ~40 min smoke run (8 steps of
300 s on a 32 x 12 x 32 shell; the latitude count is pinned near 32 so the
two free-slip ghost rows past +-80 deg stay below the pole, J = r^2 cos phi >
0):

* nothing goes non-finite and the zonal jet stays bounded (no blow-up);
* the zonal-wind field barely drifts (l2 change << the jet amplitude, and
  << the eventual O(10 m/s) baroclinic-wave perturbation);
* the case generates essentially no meridional wind — the balance tripwire.
  Measured: the jet holds at ~27.5 m/s (drift < 0.02 m/s), the geostrophic
  adjustment saturates the meridional wind at ~0.18 m/s (0.7% of the jet) and
  then plateaus. A wrong Coriolis sign/factor (the mirrored-embedding
  pseudovector trap) or a mis-placed meridional pressure gradient would drive
  O(jet) meridional wind within a few steps, as it does for TC2.

Coarse numpy smoke; the production multi-day device run is pt 3.
"""

import logging

import numpy as np
import pytest

logging.disable(logging.INFO)

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import test_hj_baroclinic as hj
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.io.debug import NullDebugWriter


def _build(nx, ny, nz, nsteps, dt=None, proj=False):
    udo = hj.UserData()
    udo.inx, udo.iny, udo.inz = nx + 1, ny + 1, nz + 1
    udo.stepmax = nsteps
    udo.initial_projection = proj
    if dt is not None:
        udo.dtfixed = udo.dtfixed0 = dt
    udo.diag = False
    udo.output_timesteps = False
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = hj.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def _winds_ms(mem, ud):
    """Interior zonal (e_lambda) and meridional (e_phi) wind fields [m/s]."""
    e = mem.elem
    i2 = (slice(2, -2), slice(2, -2), slice(2, -2))
    lam = e.x[2:-2].reshape(-1, 1, 1)
    phi = e.z[2:-2].reshape(1, 1, -1)
    cl, sl = np.cos(lam), np.sin(lam)
    cp, sp = np.cos(phi), np.sin(phi)
    e_lam = (-sl, cl, 0.0 * lam + 0.0 * phi)
    e_phi = (-sp * cl, -sp * sl, -cp + 0.0 * lam)
    rho = mem.sol.rho[i2]
    ru, rv, rw = mem.sol.rhou[i2], mem.sol.rhov[i2], mem.sol.rhow[i2]
    proj = lambda ev: (ev[0] * ru + ev[1] * rv + ev[2] * rw) / rho
    return proj(e_lam) * ud.u_ref, proj(e_phi) * ud.u_ref


def _run(mem, ud):
    mem = time_update.do(
        mem, ud, tout=1.0e9, bld=None, writer=None, debug_writer=NullDebugWriter()
    )
    assert mem.time.step == ud.stepmax
    return mem


def test_background_stays_steady():
    """No topography -> the balanced Ullrich background barely moves."""
    mem, ud = _build(32, 12, 32, nsteps=8)
    u0, v0 = _winds_ms(mem, ud)
    jet = np.abs(u0).max()
    # the analytic momenta are purely zonal: no meridional wind at init
    assert np.abs(v0).max() < 1e-10, np.abs(v0).max()
    assert 25.0 < jet < 30.0, jet  # ~28 m/s Ullrich jet (Fig. 1a)

    mem = _run(mem, ud)
    u1, v1 = _winds_ms(mem, ud)

    for a in ("rho", "rhou", "rhov", "rhow", "rhoY"):
        assert np.all(np.isfinite(getattr(mem.sol, a))), a

    zonal_drift = np.sqrt(np.mean((u1 - u0) ** 2))
    merid_wind = np.abs(v1).max()

    # jet stays bounded (no blow-up): it holds to ~0.05% (measured drift <0.02)
    assert u1.max() < 1.05 * jet, u1.max()
    # zonal field barely drifts, and stays << the eventual wave amplitude
    # (measured ~0.012 m/s = 0.04% of the jet; gate at 25x headroom)
    assert zonal_drift < 0.01 * jet, zonal_drift
    # balance tripwire: the geostrophic adjustment saturates at ~0.18 m/s
    # (0.7% of the jet); a broken balance would drive O(jet). Gate at 0.05*jet
    # (~1.4 m/s) -> ~8x over the adjustment, ~20x under a real imbalance.
    assert merid_wind < 0.05 * jet, merid_wind
