"""Terrain h == 0 identity oracle.

Two guarantees, two mechanisms:

1. **Bypass identity (bit-exact).** No registered golden-master case
   defines ``orography``, so ``elem.metric is None`` and every solver
   call site takes the uniform-Cartesian path untouched. The 9-case
   golden-master suite (``test_flow_solver.py``) is the authority; this
   file just pins the precondition.

2. **Forced-flat identity (~1e-13).** The metric machinery switched ON
   with ``h == 0`` is algebraically the identity (J == 1, G == 0) but
   multiplies extra factors through the operators. As each metric-aware
   operator lands, a comparison of the flat-metric path against the
   uniform-Cartesian path is added here — run in both the compressible
   and pseudo-incompressible regimes to cover the pi-update / wcenter
   consistency factors.
"""

import importlib

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.interfaces.ic_config import IC_MODULES
from pybella.utils import user_data
from pybella.utils import options as opts
from pybella.utils.operators import divergence

# the golden-master regression cases (smoke_* cases are exempt: they own
# no target and are allowed to grow terrain)
GOLDEN_MASTER_CASES = [
    "test_travelling_vortex",
    "test_travelling_vortex_3d_coriolis",
    "test_internal_long_wave",
    "test_igw_baldauf_brdar",
    "test_lamb_wave",
    "test_blending_warm_bubble",
    "test_unstable_lamb",
    "test_swe_vortex",
    "test_straka",
]


def test_golden_master_cases_have_no_orography():
    """Precondition for bit-identity: every target-bearing case must take
    the metric-bypass path (elem.metric is None)."""
    for ic in GOLDEN_MASTER_CASES:
        module = importlib.import_module(IC_MODULES[ic])
        ud = user_data.UserDataInit(**vars(module.UserData()))
        assert getattr(ud, "orography", None) is None, (
            f"{ic} defines orography — golden-master cases must stay on "
            "the uniform-Cartesian bypass path"
        )


# --- forced-flat operator comparisons (one per metric-aware operator) ------
#
# Not yet covered here: elliptic operator C_ij and wcenter; advective
# fluxes + CFL.


class _StubUD:
    def __init__(self, v=1, orography=None):
        self.inx, self.iny, self.inz = 13, 9, 7
        self.xmin, self.xmax = -1.0, 1.0
        self.ymin, self.ymax = 0.0, 2.0
        self.zmin, self.zmax = -0.5, 0.5
        self.bdry_type = np.array(
            [opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.PERIODIC]
        )
        self.gravity_direction = v
        if orography is not None:
            self.orography = orography


class _StubSol:
    """Smooth, deterministic 3D fields (no randomness, reproducible)."""

    def __init__(self, elem):
        # x + y + z forces every field to the full (icx, icy, icz) shape
        x, y, z = np.meshgrid(elem.x, elem.y, elem.z, indexing="ij")
        self.rho = 1.0 + 0.1 * np.sin(x) * np.cos(y) + 0.05 * z
        self.rhoY = self.rho * (1.0 + 0.02 * np.cos(x + y - z))
        self.rhou = np.sin(2 * x) + 0.3 * y * z
        self.rhov = np.cos(y) * (1.0 + 0.2 * x) + 0.0 * z
        self.rhow = 0.5 * np.sin(z + x) + 0.0 * y


def _rhs_for(ud_metric):
    elem, node = dis_grid.grid_init(ud_metric)
    sol = _StubSol(elem)
    rhs = np.zeros(node.isc)
    divergence.compute_at_nodes(rhs, elem, sol, ud_metric)
    return rhs


def test_divergence_flat_metric_matches_plain():
    flat = lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2
    rhs_plain = _rhs_for(_StubUD())
    rhs_flat = _rhs_for(_StubUD(orography=flat))
    assert np.max(np.abs(rhs_flat - rhs_plain)) <= 1e-13


def test_gradient_map_flat_metric_is_exact_identity():
    flat = lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2
    elem, _ = dis_grid.grid_init(_StubUD(orography=flat))
    from pybella.flow_solver.discretisation import terrain

    sol = _StubSol(elem)
    dp0 = [sol.rhou.copy(), sol.rhov.copy(), sol.rhow.copy()]
    dp = terrain.apply_gradient_map(elem.metric, [d.copy() for d in dp0])
    for a, b in zip(dp, dp0):
        # J == 1, G == 0: multiplication by 1.0 and subtraction of 0.0 are
        # exact in floating point — bit-identity, not approximation
        assert np.array_equal(a, b)


@pytest.mark.parametrize("v", [0, 1, 2])
def test_gradient_map_constant_metric_algebra(v):
    from pybella.utils import axes
    from pybella.flow_solver.discretisation.terrain import (
        MetricFields,
        apply_gradient_map,
    )

    elem, _ = dis_grid.grid_init(_StubUD(v=v))
    sol = _StubSol(elem)
    dp0 = [sol.rhou.copy(), sol.rhov.copy(), sol.rhow.copy()]
    a_h1, a_h2 = axes.horizontal_axes(v)
    shape = dp0[0].shape
    J0, g1, g2 = 2.0, 0.3, -0.7
    metric = MetricFields(
        J=np.full(shape, J0),
        G1=np.full(shape, g1),
        G2=np.full(shape, g2),
        z=np.zeros(shape),
        vaxis=v,
        haxes=(a_h1, a_h2),
    )
    dp = apply_gradient_map(metric, [d.copy() for d in dp0])
    assert np.allclose(dp[a_h1], dp0[a_h1] - (g1 / J0) * dp0[v], atol=1e-15)
    assert np.allclose(dp[a_h2], dp0[a_h2] - (g2 / J0) * dp0[v], atol=1e-15)
    assert np.allclose(dp[v], dp0[v] / J0, atol=1e-15)


@pytest.mark.parametrize("v", [0, 1, 2])
def test_divergence_metric_wiring(v):
    """Algebraic wiring checks of the role-axis mapping, exact to roundoff.

    With theta == 1 (rho == rhoY), constant J and constant slopes, the
    metric divergence must equal the plain divergence of analytically
    pre-transformed momenta.
    """
    from pybella.utils import axes

    ud = _StubUD(v=v)
    ud.bdry_type = np.array([opts.BdryType.PERIODIC] * 3)
    elem, node = dis_grid.grid_init(ud)
    sol = _StubSol(elem)
    sol.rhoY = sol.rho.copy()  # theta == 1

    a_h1, a_h2 = axes.horizontal_axes(v)
    moms = (sol.rhou, sol.rhov, sol.rhow)

    # hand-build a constant metric (bypasses the transform on purpose)
    from pybella.flow_solver.discretisation.terrain import MetricFields

    shape = sol.rho.shape
    J0, g1, g2 = 2.0, 0.3, -0.7
    elem.metric = MetricFields(
        J=np.full(shape, J0),
        G1=np.full(shape, g1),
        G2=np.full(shape, g2),
        z=np.zeros(shape),
        vaxis=v,
        haxes=(a_h1, a_h2),
    )
    rhs_metric = np.zeros(node.isc)
    divergence.compute_at_nodes(rhs_metric, elem, sol, ud)
    elem.metric = None

    # reference: plain divergence of the pre-transformed momenta
    ref_moms = [None, None, None]
    ref_moms[a_h1] = J0 * moms[a_h1]
    ref_moms[a_h2] = J0 * moms[a_h2]
    ref_moms[v] = moms[v] - g1 * moms[a_h1] - g2 * moms[a_h2]
    sol.rhou, sol.rhov, sol.rhow = ref_moms
    rhs_ref = np.zeros(node.isc)
    divergence.compute_at_nodes(rhs_ref, elem, sol, ud)

    assert np.max(np.abs(rhs_metric - rhs_ref)) <= 1e-13


# --- native 2D (lap2D path) -------------------------------------------------


class _StubSol2D:
    """Smooth, deterministic 2D fields."""

    def __init__(self, elem):
        x, y = np.meshgrid(elem.x, elem.y, indexing="ij")
        self.rho = 1.0 + 0.1 * np.sin(x) * np.cos(y)
        self.rhoY = self.rho * (1.0 + 0.02 * np.cos(x + y))
        self.rhou = np.sin(2 * x) + 0.3 * y
        self.rhov = np.cos(y) * (1.0 + 0.2 * x)
        self.rhow = np.zeros_like(x)


def _rhs_for_2d(ud_metric):
    ud_metric.inz = 1
    elem, node = dis_grid.grid_init(ud_metric)
    sol = _StubSol2D(elem)
    rhs = np.zeros((elem.icx - 1, elem.icy - 1))
    divergence.compute_at_nodes(rhs, elem, sol, ud_metric)
    return rhs


def test_divergence_2d_flat_metric_matches_plain():
    flat = lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2
    rhs_plain = _rhs_for_2d(_StubUD())
    rhs_flat = _rhs_for_2d(_StubUD(orography=flat))
    assert np.max(np.abs(rhs_flat - rhs_plain)) <= 1e-13


def test_divergence_2d_metric_wiring():
    """2D sibling of test_divergence_metric_wiring (haxes = (0, None))."""
    from pybella.flow_solver.discretisation.terrain import MetricFields

    ud = _StubUD()
    ud.inz = 1
    ud.bdry_type = np.array([opts.BdryType.PERIODIC] * 3)
    elem, node = dis_grid.grid_init(ud)
    sol = _StubSol2D(elem)
    sol.rhoY = sol.rho.copy()  # theta == 1

    shape = sol.rho.shape
    J0, g1 = 2.0, 0.3
    elem.metric = MetricFields(
        J=np.full(shape, J0),
        G1=np.full(shape, g1),
        G2=None,
        z=np.zeros(shape),
        vaxis=1,
        haxes=(0, None),
    )
    rhs_metric = np.zeros((elem.icx - 1, elem.icy - 1))
    divergence.compute_at_nodes(rhs_metric, elem, sol, ud)
    elem.metric = None

    rhou0 = sol.rhou.copy()
    sol.rhou = J0 * rhou0
    sol.rhov = sol.rhov - g1 * rhou0
    rhs_ref = np.zeros((elem.icx - 1, elem.icy - 1))
    divergence.compute_at_nodes(rhs_ref, elem, sol, ud)

    assert np.max(np.abs(rhs_metric - rhs_ref)) <= 1e-13
