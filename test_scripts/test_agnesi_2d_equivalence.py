"""Native-2D vs quasi-2D-3D Agnesi equivalence — the lap2D-terrain payoff proof.

Two gates:

1. **Path comparison** (15 golden-master steps, 128x64): native 2D (lap2D
   with terrain-folded cross terms) vs the proven quasi-2D 3D path (lap3D
   full tensor). The two laplacian families carry *historically different
   wall discretizations* — measured on a FLAT wall-bounded impulse, the
   pre-existing 2D-vs-3D gap is already ~13% rel-L2 in w / ~19% in p2
   after 5 steps (2026-06-10, solver-tolerance independent). The gates
   here sit just above that honest cross-discretization floor: they catch
   sign/axis/J-factor errors (which blow up by orders of magnitude), not
   stencil-convention differences the flat solver already had.

2. **Absolute physics** (the strong gate): the native-2D run passes the
   same Smith (1980) analytic-oracle thresholds as the 3D path —
   wrongness, not just change. Same reduced config as
   ``test_agnesi_analytic`` (96x48, t U / a = 12), roughly halved runtime
   in 2D. Calibration 2026-06-10 (native 2D): w 0.40, u' 0.43,
   drag_ratio 0.98, flux_constancy 0.04 — within a percent of the 3D
   path's values, which is the real equivalence statement.
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


def _run(inz, steps=None, grid=None):
    ud = user_data.UserDataInit(**vars(case.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.inz = inz
    ud.tout = [1e6]
    ud.diag = False
    if steps is not None:
        ud.stepmax = steps
    if grid is not None:
        ud.inx, ud.iny, ud.inbcy = grid
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = case.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter()), ud


def _slab(arr, ndim):
    if ndim == 2:
        return arr[2:-2, 2:-2]
    return arr[2:-2, 2:-2, 0]


def test_native_2d_matches_quasi_2d_3d():
    mem3, _ = _run(inz=2)
    mem2, _ = _run(inz=1)

    # measured 2026-06-10: rho 4.9e-6, rhou 9.4e-5, rhov 8.8e-2 (wave
    # perturbation on a near-zero field), rhoY 4.9e-6, rhoX 1.1e-3
    for attr, tol in (
        ("rho", 2e-5),
        ("rhou", 4e-4),
        ("rhov", 0.2),
        ("rhoY", 2e-5),
        ("rhoX", 5e-3),
    ):
        a = _slab(getattr(mem3.sol, attr), mem3.elem.ndim)
        b = _slab(getattr(mem2.sol, attr), mem2.elem.ndim)
        scale = max(np.linalg.norm(a), 1e-30)
        err = np.linalg.norm(a - b) / scale
        assert err <= tol, f"{attr}: 2D vs 3D rel-L2 {err:.3e} > {tol}"

    # wave field: measured 0.11 (the wall-convention floor); an axis/sign/J
    # defect lands at O(1)
    w3 = _slab(mem3.sol.rhov / mem3.sol.rho, mem3.elem.ndim)
    w2 = _slab(mem2.sol.rhov / mem2.sol.rho, mem2.elem.ndim)
    err_w = np.linalg.norm(w3 - w2) / np.linalg.norm(w3)
    assert err_w <= 0.25, f"w wave field: 2D vs 3D rel-L2 {err_w:.3e}"


def test_native_2d_passes_smith_oracle():
    """Native 2D vs the Smith (1980) analytic solution — same gates as the
    3D oracle in test_agnesi_analytic.py."""
    mem, ud = _run(inz=1, steps=240, grid=(96 + 1, 48 + 1, 24))
    metrics, _ = smith.compare(mem, ud, z_lo_SI=1000.0, z_hi_SI=4500.0)

    assert metrics["w"] <= 0.55, f"w rel-L2 vs Smith: {metrics['w']:.3f}"
    assert metrics["u"] <= 0.60, f"u' rel-L2 vs Smith: {metrics['u']:.3f}"
    assert (
        0.85 <= metrics["drag_ratio"] <= 1.15
    ), f"momentum flux / analytic drag: {metrics['drag_ratio']:.3f}"
    assert (
        metrics["flux_constancy"] <= 0.10
    ), f"momentum-flux height variation: {metrics['flux_constancy']:.3f}"
