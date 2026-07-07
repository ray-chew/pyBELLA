"""Stage F, increment F6: TC6 Rossby-Haurwitz on the full pole-to-pole sphere.

TC6 has real longitude structure at all latitudes, so this exercises the
pole ghost exchange (F1), pole-face flux closure (F2) and polar filter (F3)
under genuine dynamics (TC2 is longitude-independent). Short tripwires (the
7-day RH-4 phase-speed validation runs from run_scripts):

- the wave stays bounded and finite over 10 steps (no pole blow-up);
- J-weighted mass is conserved to machine precision (the pole-face closure
  under full dynamics);
- the pattern keeps its RH-4 structure (dominant zonal wavenumber 4).
"""

import logging

import numpy as np

logging.disable(logging.INFO)

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.tests import test_sphere_swe_tc6_global as tc6g
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.io.debug import NullDebugWriter


def _run(nsteps, nlam=48, nphi=54):
    udo = tc6g.UserData()
    udo.inx, udo.inz = nlam + 1, nphi + 1
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
    sol = tc6g.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())

    inner = (slice(2, -2), slice(2, -2), slice(2, -2))
    J = elem.metric.J
    mass0 = np.sum((J * sol.rho)[inner])
    h0 = sol.rho[inner].copy()
    mem = time_update.do(
        mem, ud, tout=1.0e9, bld=None, writer=None, debug_writer=NullDebugWriter()
    )
    assert mem.time.step == nsteps
    mass1 = np.sum((J * mem.sol.rho)[inner])
    return mem, inner, h0, mass0, mass1


def test_tc6_global_bounded_conservative_rh4():
    mem, inner, h0, mass0, mass1 = _run(10)
    h = mem.sol.rho[inner]
    assert np.isfinite(mem.sol.rho).all()

    # bounded evolution (RH-4 wave, not steady): depth stays near its
    # initial range, no pole blow-up
    assert h.min() > 0.5 * h0.min()
    assert h.max() < 1.5 * h0.max()

    # J-weighted mass conserved to machine precision (pole-face closure)
    assert abs(mass1 - mass0) / abs(mass0) < 1e-12, (mass0, mass1)

    # the pattern keeps its RH-4 signature: the dominant zonal wavenumber of
    # the depth anomaly at mid-latitude is 4
    igl = 2
    mid = h.shape[2] // 2  # equatorial-ish latitude band
    row = mem.sol.rho[igl:-igl, mem.sol.rho.shape[1] // 2, mid + inner[2].start]
    amp = np.abs(np.fft.rfft(row - row.mean()))
    assert np.argmax(amp) == 4, np.argmax(amp)
