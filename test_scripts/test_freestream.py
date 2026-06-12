"""Freestream-preservation tripwire for the general curvilinear metric.

A uniform flow (constant rho, rhoY, velocity) has zero physical
divergence, so the metric divergence rhs = Sum_a D_a (N_a . f) with
f = theta * m = const must vanish up to the discrete metric defect
c . Sum_a D_a N_a — the discrete analogue of the metric identity
Sum_a d_a N_a = 0 (Klein's p. 9 compatibility question). With analytic
collocated metrics (decision 1, option 1) the identity holds to
truncation only; this test measures the defect directly on a wavy,
genuinely 3D Tier-2 map and trips if it stops converging at second
order or exceeds an absolute bound.

If a future strongly distorted (Tier-3) grid trips this, the planned
fallback is discrete face-based normals (freestream-exact by
construction) — see tfc_generalization_plan.md, Phase 2.
"""

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import terrain
from pybella.flow_solver.utils import fields
from pybella.utils import options as opts
from pybella.utils.operators import divergence

pytestmark = pytest.mark.skipif(
    not hasattr(terrain, "build_metric_fields_from_map"),
    reason="general metric machinery (tfc pt 1) not present",
)


class _StubUD:
    """Minimal ud: all-PERIODIC so compute_at_nodes zeroes no wall slabs
    (the uniform flow must reach the stencil unclipped)."""

    def __init__(self, n):
        self.inx, self.iny, self.inz = n, n, n
        self.xmin, self.xmax = -1.0, 1.0
        self.ymin, self.ymax = 0.0, 2.0
        self.zmin, self.zmax = -0.8, 0.8
        self.bdry_type = np.array([opts.BdryType.PERIODIC] * 3)
        self.gravity_direction = 1
        # terrain must be "active" for grid_init; the test overrides the
        # metric with the wavy general map below
        self.orography = lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2


class _WavyMap(terrain.CurvilinearMap):
    """Genuinely 3D Tier-2 map: stretched x and y, wavy Gal-Chen z."""

    eta0, etat = 0.0, 2.0

    def _h(self, xi1, xi2):
        return 0.08 * np.cos(np.pi * xi1) * (1.0 + 0.4 * np.sin(2.5 * xi2))

    def _h1(self, xi1, xi2):
        return -0.08 * np.pi * np.sin(np.pi * xi1) * (1.0 + 0.4 * np.sin(2.5 * xi2))

    def _h2(self, xi1, xi2):
        return 0.08 * np.cos(np.pi * xi1) * 0.4 * 2.5 * np.cos(2.5 * xi2)

    def _decay(self, eta):
        return (self.etat - eta) / (self.etat - self.eta0)

    def coordinates(self, xi):
        xi1, eta, xi2 = xi
        z = eta + self._h(xi1, xi2) * self._decay(eta)
        x = xi1 + 0.1 * np.sin(np.pi * xi1)
        y = xi2 + 0.08 * np.sin(1.25 * np.pi * xi2)
        return [x + 0.0 * eta + 0.0 * xi2, z, y + 0.0 * xi1 + 0.0 * eta]

    def tangents(self, xi):
        xi1, eta, xi2 = xi
        b = self._decay(eta)
        J = (1.0 - self._h(xi1, xi2) / (self.etat - self.eta0)) + 0.0 * eta
        G1 = self._h1(xi1, xi2) * b
        G2 = self._h2(xi1, xi2) * b
        xp = 1.0 + 0.1 * np.pi * np.cos(np.pi * xi1)
        yp = 1.0 + 0.08 * 1.25 * np.pi * np.cos(1.25 * np.pi * xi2)
        zero = 0.0 * (J + xp + yp)
        one = 1.0 + zero
        return [
            [xp + zero, G1 + zero, zero],
            [zero, J + zero, zero],
            [zero, G2 + zero, yp + zero],
        ]


def _freestream_defect(n):
    """Max relative divergence of a uniform flow on the wavy n^3 grid."""
    ud = _StubUD(n)
    elem, node = dis_grid.grid_init(ud)
    elem.metric = terrain.build_metric_fields_from_map(elem, ud, _WavyMap())

    sol = fields.CellSolField(elem.sc)
    sol.rho[...] = 1.0
    sol.rhoY[...] = 1.2
    u0, v0, w0 = 0.7, -0.4, 0.3
    sol.rhou[...] = u0
    sol.rhov[...] = v0
    sol.rhow[...] = w0

    rhs = np.zeros(node.isc)
    divergence.compute_at_nodes(rhs, elem, sol, ud)

    # normalize by the flux scale over the smallest spacing
    theta = 1.2
    fmax = (
        theta
        * max(abs(u0), abs(v0), abs(w0))
        * max(np.max(np.abs(c)) for Na in elem.metric.N for c in Na)
    )
    return np.max(np.abs(rhs)) / (fmax / min(elem.dxyz))


def test_freestream_defect_small_and_second_order():
    d_coarse = _freestream_defect(16)
    d_fine = _freestream_defect(32)
    # absolute tripwire: the defect is a truncation error, not O(1)
    assert d_coarse < 5e-3, f"freestream defect {d_coarse:.2e} too large"
    # convergence tripwire: second-order metrics must shrink ~4x per halving
    assert (
        d_fine < d_coarse / 3.0
    ), f"freestream defect not converging: {d_coarse:.2e} -> {d_fine:.2e}"


def test_freestream_exact_on_flat_metric():
    """h == 0, no stretch: N_a = e_a exactly, so a uniform flow yields a
    bit-exact zero divergence (the bypass analogue on the general path)."""
    ud = _StubUD(12)
    elem, node = dis_grid.grid_init(ud)
    # ud.orography is flat -> builder gives the identity metric
    assert np.all(elem.metric.J == 1.0)

    sol = fields.CellSolField(elem.sc)
    sol.rho[...] = 1.0
    sol.rhoY[...] = 1.0
    sol.rhou[...] = 0.5
    sol.rhov[...] = 0.5
    sol.rhow[...] = 0.5

    rhs = np.zeros(node.isc)
    divergence.compute_at_nodes(rhs, elem, sol, ud)
    assert np.all(rhs == 0.0)
