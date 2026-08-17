"""Resting atmosphere over terrain — the classic TFC physics gate.

A hydrostatically balanced isothermal atmosphere at rest over a witch-of-
Agnesi hill must stay (nearly) at rest: every spurious metric force —
imbalanced pressure-gradient mapping, wrong Jacobian factors, ghost-cell
hydrostatics at the wrong height — shows up as wind generated from nothing.
Discrete rest is not exact (corner-averaged gradients across columns of
different physical height carry truncation error), so the gate is a small
velocity bound, not machine zero.

Also pins the h == 0 equivalence of the field-mode hydrostates against the
1D column profiles (exact for the analytical state, quadrature-tight for
the integrated state).
"""

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import hydrostatics, thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import smoke_agnesi
from pybella.utils import axes, user_data
from pybella.utils.data_structures import ModelState


class _StubWriter:
    def write(self, *args, **kwargs):
        pass

    def populate(self, *args, **kwargs):
        pass

    def write_all(self, *args, **kwargs):
        pass


def _agnesi(ud, h0_m=400.0, a_m=5000.0):
    h0, a = h0_m / ud.h_ref, a_m / ud.h_ref
    return lambda xi1, xi2: h0 * a**2 / (xi1**2 + a**2) + 0.0 * xi2


def _make_resting_mem(orography=True, steps=5, inz=2, sleve=False):
    ud = user_data.UserDataInit(**vars(smoke_agnesi.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.u_wind_speed = 0.0
    ud.stepmax = steps
    ud.inz = inz  # 2 = quasi-2D 3D (lap3D path), 1 = native 2D (lap2D path)
    ud.tout = [1e6]  # step-limited
    # smoke_agnesi carries its own hill; override per variant
    ud.orography = _agnesi(ud) if orography else None
    if sleve:
        from pybella.flow_solver.discretisation import terrain

        hill = ud.orography
        ud.orography_smooth = lambda xi1, xi2: 0.5 * hill(xi1, xi2)
        ud.vertical_transform = terrain.SLEVETransform(
            s1=6000.0 / ud.h_ref, s2=1500.0 / ud.h_ref
        )
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = smoke_agnesi.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def _max_speed_ms(mem, ud):
    i2 = tuple(slice(2, -2) for _ in range(mem.elem.ndim))
    speeds = [
        np.max(np.abs(getattr(mem.sol, m)[i2] / mem.sol.rho[i2]))
        for m in ("rhou", "rhov", "rhow")
    ]
    return max(speeds) * ud.u_ref


@pytest.mark.parametrize("inz", [2, 1], ids=["q2d3d", "native2d"])
def test_resting_atmosphere_over_hill(inz):
    mem, ud = _make_resting_mem(orography=True, steps=5, inz=inz)
    mem = time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())
    vmax = _max_speed_ms(mem, ud)
    # measured ~1e-10 m/s (the bicgstab solve floor, same as flat): with the
    # p2 == 0 perturbation-pressure convention the discrete rest state is
    # exact through the metric machinery, not merely truncation-small
    assert vmax < 1e-8, f"spurious wind over terrain: max |v| = {vmax:.3e} m/s"


@pytest.mark.parametrize("inz", [2, 1], ids=["q2d3d", "native2d"])
def test_resting_atmosphere_over_hill_sleve(inz):
    """Balanced rest under the first eta-dependent Jacobian: field-mode
    hydrostates at SLEVE heights, bottom BC with J varying through the
    ghost rows, metric advection and elliptic assembly all see J(eta)."""
    mem, ud = _make_resting_mem(orography=True, steps=5, inz=inz, sleve=True)
    mem = time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())
    vmax = _max_speed_ms(mem, ud)
    assert vmax < 1e-8, f"spurious wind over SLEVE terrain: {vmax:.3e} m/s"


@pytest.mark.parametrize("inz", [2, 1], ids=["q2d3d", "native2d"])
def test_resting_atmosphere_flat_stays_still(inz):
    """Flat baseline: quantifies the no-terrain spurious-wind floor."""
    mem, ud = _make_resting_mem(orography=False, steps=5, inz=inz)
    mem = time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())
    vmax = _max_speed_ms(mem, ud)
    assert vmax < 1e-8, f"flat resting atmosphere drifted: {vmax:.3e} m/s"


@pytest.mark.parametrize("inz", [2, 1], ids=["q2d3d", "native2d"])
def test_uniform_flow_flat_metric_matches_plain(inz):
    """Full time loop, wind on: forced-flat metric == plain to roundoff.

    Exercises every metric-aware branch (divergence, gradients, elliptic,
    ghost cells, sweep-oriented metric flips) with J == 1, G == 0 against
    the untouched uniform-Cartesian path.
    """

    def run(orography):
        ud = user_data.UserDataInit(**vars(smoke_agnesi.UserData()))
        ud.coriolis_strength = np.array(ud.coriolis_strength)
        ud.stepmax = 3
        ud.tout = [1e6]
        ud.inz = inz
        # forced-flat metric vs true bypass (the case's own hill is overridden)
        ud.orography = (lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2) if orography else None
        elem, node = dis_grid.grid_init(ud)
        sol = fields.CellSolField(elem.sc)
        th = thermodynamics.ThermodynamicalQuantities(ud)
        npf = fields.NodePressureField(elem, node, ud)
        sol = smoke_agnesi.sol_init(sol, npf, elem, node, th, ud)
        mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
        bdry_c.set_ghost_cells(mem, ud)
        return time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())

    plain = run(False)
    flat = run(True)
    for attr in ("rho", "rhou", "rhov", "rhow", "rhoY"):
        a = getattr(plain.sol, attr)
        b = getattr(flat.sol, attr)
        scale = max(np.max(np.abs(a)), 1.0)
        err = np.max(np.abs(a - b)) / scale
        assert err <= 1e-12, f"{attr}: rel {err:.2e}"


@pytest.mark.parametrize("inz", [2, 1], ids=["q2d3d", "native2d"])
def test_mountain_wave_smoke_conservation_and_response(inz):
    """First end-to-end terrain run: 10 m/s wind over the 400 m hill.

    Gates: (i) J-weighted mass and P = rho*Y are conserved by the metric
    advection (flux form, periodic x, no-flux walls), (ii) the hill
    actually forces a vertical-velocity response of a sane magnitude
    (linear estimate U * max|dh/dx| ~ 0.5 m/s for these parameters).
    """
    ud = user_data.UserDataInit(**vars(smoke_agnesi.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.stepmax = 5
    ud.tout = [1e6]
    ud.inz = inz
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = smoke_agnesi.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)

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
    assert 1e-3 < w_ms < 5.0, f"mountain-wave w response out of range: {w_ms:.3e} m/s"


def test_field_mode_hydrostates_match_profiles_when_flat():
    flat = lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2

    def states_pair(builder):
        ud0 = user_data.UserDataInit(**vars(smoke_agnesi.UserData()))
        ud0.coriolis_strength = np.array(ud0.coriolis_strength)
        ud0.orography = None
        elem0, node0 = dis_grid.grid_init(ud0)
        npf0 = fields.NodePressureField(elem0, node0, ud0)
        builder(npf0, elem0, node0, thermodynamics.ThermodynamicalQuantities(ud0), ud0)

        ud1 = user_data.UserDataInit(**vars(smoke_agnesi.UserData()))
        ud1.coriolis_strength = np.array(ud1.coriolis_strength)
        ud1.orography = flat
        elem1, node1 = dis_grid.grid_init(ud1)
        npf1 = fields.NodePressureField(elem1, node1, ud1)
        builder(npf1, elem1, node1, thermodynamics.ThermodynamicalQuantities(ud1), ud1)

        v = axes.vertical_axis(ud0)
        return npf0, npf1, elem0, v

    # analytical_state: identical expressions either way -> exact.
    # integrated_state: the 1D profile uses coarse per-cell trapezoids,
    # the field branch a fine-grid quadrature — they differ by the COARSE
    # scheme's truncation (~1.5e-5 here); the bound guards wiring errors
    # (axis, sign, reference level), not quadrature equivalence.
    for builder, tol in (
        (hydrostatics.analytical_state, 0.0),
        (hydrostatics.integrated_state, 5e-5),
    ):
        npf0, npf1, elem0, v = states_pair(builder)
        for attr in ("rho0", "rhoY0", "p0", "p20", "S0", "Y0"):
            profile = axes.expand_profile(
                getattr(npf0.HydroState, attr), elem0.ndim, v, elem0.sc
            )
            field = getattr(npf1.HydroState, attr)
            err = np.max(np.abs(field - profile)) / np.max(np.abs(profile))
            assert err <= tol, f"{builder.__name__}.{attr}: rel {err:.2e} > {tol}"
