"""Tier-2 physics gates: terrain + horizontal stretching, end to end.

First full runs on a genuinely curvilinear (general-path) grid: the
x-stretched Gal-Chen map exercises the N-flux advection, the A-mapped
gradients, the general elliptic fold, the contravariant CFL, the
effective-slope wall reflection and the z_eta-based ghost/hydrostate
thickness all at once.

Gates (the vertical-line smoke bars, re-pointed at Tier 2):
1. resting atmosphere stays at rest (< 1e-8 m/s — the p2 == 0
   perturbation convention makes discrete rest metric-independent),
2. mountain-wave smoke: J-weighted mass and P = rho*Y conserved to
   1e-10, vertical response of sane magnitude.
"""

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import terrain
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import smoke_agnesi
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState

pytestmark = pytest.mark.skipif(
    not hasattr(terrain, "build_metric_fields_from_map"),
    reason="general metric machinery (tfc pt 1) not present",
)


class _StubWriter:
    def write(self, *args, **kwargs):
        pass

    def populate(self, *args, **kwargs):
        pass

    def write_all(self, *args, **kwargs):
        pass


class _StretchedHillMap(terrain.CurvilinearMap):
    """Periodic x-stretch + periodic hill, Gal-Chen z (Tier-2 map).

    Periodic-consistent in xi1: stretch displacement and hill share the
    domain period, so the metric is smooth across the x seam.
    """

    def __init__(self, ud, h0_m=400.0, ax_rel=0.05):
        self.L = ud.xmax - ud.xmin
        self.h0 = h0_m / ud.h_ref
        self.ax = ax_rel * self.L
        self.eta0, self.etat = ud.ymin, ud.ymax

    def _h(self, xi1):
        return self.h0 * np.cos(np.pi * xi1 / self.L) ** 2

    def _dh(self, xi1):
        return -self.h0 * (np.pi / self.L) * np.sin(2.0 * np.pi * xi1 / self.L)

    def _decay(self, eta):
        return (self.etat - eta) / (self.etat - self.eta0)

    def coordinates(self, xi):
        xi1, eta = xi[0], xi[1]
        x = xi1 + self.ax * np.sin(2.0 * np.pi * xi1 / self.L)
        z = eta + self._h(xi1) * self._decay(eta)
        if len(xi) == 2:
            return [x + 0.0 * eta, z]
        return [x + 0.0 * eta + 0.0 * xi[2], z, None]

    def tangents(self, xi):
        xi1, eta = xi[0], xi[1]
        xp = 1.0 + self.ax * (2.0 * np.pi / self.L) * np.cos(2.0 * np.pi * xi1 / self.L)
        J = (1.0 - self._h(xi1) / (self.etat - self.eta0)) + 0.0 * eta
        G1 = self._dh(xi1) * self._decay(eta)
        zero = 0.0 * (J + xp)
        if len(xi) == 2:
            return [[xp + zero, G1 + zero], [zero, J + zero]]
        zero3 = zero + 0.0 * xi[2]
        one = 1.0 + zero3
        return [
            [xp + zero3, G1 + zero3, zero3],
            [zero3, J + zero3, zero3],
            [zero3, zero3, one],
        ]


def _make_mem(wind=True, steps=5, inz=2):
    ud = user_data.UserDataInit(**vars(smoke_agnesi.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    if not wind:
        ud.u_wind_speed = 0.0
    ud.stepmax = steps
    ud.inz = inz
    ud.tout = [1e6]  # step-limited
    ud.orography = None  # the general map below replaces the legacy builder
    elem, node = dis_grid.grid_init(ud)
    cmap = _StretchedHillMap(ud)
    elem.metric = terrain.build_metric_fields_from_map(elem, ud, cmap)
    node.metric = terrain.build_metric_fields_from_map(node, ud, cmap)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = smoke_agnesi.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


@pytest.mark.parametrize("inz", [2, 1], ids=["q2d3d", "native2d"])
def test_stretched_resting_atmosphere(inz):
    mem, ud = _make_mem(wind=False, steps=5, inz=inz)
    mem = time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())
    i2 = tuple(slice(2, -2) for _ in range(mem.elem.ndim))
    speeds = [
        np.max(np.abs(getattr(mem.sol, m)[i2] / mem.sol.rho[i2]))
        for m in ("rhou", "rhov", "rhow")
    ]
    vmax = max(speeds) * ud.u_ref
    assert vmax < 1e-8, f"spurious wind on stretched grid: {vmax:.3e} m/s"


@pytest.mark.parametrize("inz", [2, 1], ids=["q2d3d", "native2d"])
def test_stretched_mountain_wave_smoke(inz):
    mem, ud = _make_mem(wind=True, steps=5, inz=inz)
    elem = mem.elem

    i2 = tuple(slice(2, -2) for _ in range(elem.ndim))
    J = elem.metric.J
    mass0 = np.sum((J * mem.sol.rho)[i2])
    P0 = np.sum((J * mem.sol.rhoY)[i2])

    mem = time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())

    mass1 = np.sum((J * mem.sol.rho)[i2])
    P1 = np.sum((J * mem.sol.rhoY)[i2])
    assert abs(mass1 - mass0) / mass0 < 1e-10, f"mass drift {(mass1-mass0)/mass0:.2e}"
    assert abs(P1 - P0) / P0 < 1e-10, f"P drift {(P1-P0)/P0:.2e}"

    w_ms = np.max(np.abs(mem.sol.rhov[i2] / mem.sol.rho[i2])) * ud.u_ref
    assert 1e-3 < w_ms < 5.0, f"stretched mountain-wave response: {w_ms:.3e} m/s"
