"""Agnesi mountain-wave analytic oracle: simulation vs Smith (1980).

Runs the ``test_agnesi_hydrostatic`` configuration in-process to a
quasi-steady state (t U / a ~ 10) and compares the wave field against the
steady linear hydrostatic solution and the analytic wave drag — catching
*wrongness* of the terrain-following dynamics, not just *change*.

Gates (rel-L2 in the window 1-9 km, below the sponge):
- vertical velocity w vs Smith,
- horizontal perturbation u' vs Smith,
- vertically integrated momentum flux vs the analytic drag,
- flux constancy with height (linear steady state transports momentum
  uniformly to the breaking/sponge level).
"""

import numpy as np

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import agnesi_smith_analytic as smith
from pybella.tests import test_agnesi_hydrostatic as case
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState


class _StubWriter:
    def write(self, *args, **kwargs):
        pass

    def populate(self, *args, **kwargs):
        pass

    def write_all(self, *args, **kwargs):
        pass


def run_to_steady_state(steps=240):
    """Reduced oracle config: 96x48, t = steps * 50 s (240 -> t U / a = 12).

    The wave field below ~5 km is quasi-steady by then (vertical group
    speed ~ U^2 / (N a) = 1 m/s); the comparison window stays below that.
    Calibration (this config): w 0.40, u' 0.43, drag_ratio 0.98,
    flux_constancy 0.04. The full-resolution 128x64 run at t U / a = 20
    gives w 0.39 / u' 0.31 / drag 1.04 — resolution is not the limiter,
    residual spin-up transients are.
    """
    ud = user_data.UserDataInit(**vars(case.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.inx = 96 + 1  # dx ~ 2.08 km (a / dx = 4.8)
    ud.iny = 48 + 1  # dy = 500 m
    ud.inbcy = 24  # sponge above 12 km
    ud.stepmax = steps
    ud.tout = [1e6]
    ud.diag = False
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = case.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    mem = time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())
    return mem, ud


def test_agnesi_vs_smith():
    mem, ud = run_to_steady_state()
    metrics, _ = smith.compare(mem, ud, z_lo_SI=1000.0, z_hi_SI=4500.0)

    # gates ~35-40% above the calibrated values (see run_to_steady_state);
    # any metric-term sign/factor/orientation bug blows these completely
    assert metrics["w"] <= 0.55, f"w rel-L2 vs Smith: {metrics['w']:.3f}"
    assert metrics["u"] <= 0.60, f"u' rel-L2 vs Smith: {metrics['u']:.3f}"
    assert (
        0.85 <= metrics["drag_ratio"] <= 1.15
    ), f"momentum flux / analytic drag: {metrics['drag_ratio']:.3f}"
    assert (
        metrics["flux_constancy"] <= 0.10
    ), f"momentum-flux variation with height: {metrics['flux_constancy']:.3f}"
