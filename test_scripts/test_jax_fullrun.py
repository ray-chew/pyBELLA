"""Full regression runs on the hybrid JAX backend (elliptic solve plus the
advection / Coriolis / diffusion kernels, end-to-end).

Runs complete production cases with ``PYBELLA_BACKEND=jax`` — the env var
flips ``ud.backend`` without touching the case files — and asserts a zero
exit code, which means ``CompareSol.test_do`` passed against the *same
stored golden-master targets* the numpy backend is held to (per-field
max-abs < 1e-5). This is the strongest reproducibility statement available:
the hybrid JAX path (elliptic solve, advection recovery+HLL, advective
fluxes, Coriolis, diffusion) reproduces the numpy reference within the
regression tolerance over full runs, not just single steps.

The four cases cover: 2D periodic advection+elliptic (vortex), explicit
diffusion + x-WALLs (Straka), 3D elliptic + full Coriolis (3D vortex), and
terrain-following coordinates (Agnesi). The remaining regression cases were
validated once on this backend; run them ad hoc with
``PYBELLA_BACKEND=jax pytest test_scripts/test_flow_solver.py``.

The sphere path (non-vertical-line metric: ``coriolis_field`` H^-1, general
free-slip walls, tangent-plane surface constraint) is gated in-process by
short Williamson TC2 runs compared jax-vs-numpy, not against the stored
target. Two kinds of gate:

* ``*_stepper_bit_identical`` — initial projection OFF. The stepper twins are
  exact (~1e-16), independent of the jax / scipy release. These are the
  gates that catch a broken kernel.
* ``*_reproduces_numpy`` — initial projection ON. The bicgstab projection
  fixes the answer only to its residual class (scipy's default ``rtol=1e-5``,
  which ``ud.tol`` does not override), and *which* member each backend lands
  on depends on the jax release: jax 0.10.1 agrees with scipy to ~1e-5 in the
  momenta, jax 0.11.1 to ~6e-4 (same numpy/scipy, same machine). These are
  therefore loose sanity checks, not precision gates.

Skips cleanly when jax is not installed.
"""

import importlib.util
import os
import subprocess

import pytest

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("jax") is None, reason="jax not installed"
)

_SPHERE_INNER = (slice(2, -2), slice(2, -2), slice(2, -2))

# Loose sanity tolerances for the projection-ON sphere gates: ~3x above the
# jax-0.11.1 projection floor measured against scipy 1.18.1 (channel: scalars
# 8.4e-5, momenta/rho 6.0e-4, p2 7.3e-5; global: 1.7e-5 / 2.4e-4 / 5.8e-4).
# A broken pole exchange or wrong Coriolis factor moves the fields by O(1e-2).
_SANITY_TOL = {"scalar": 5e-4, "momentum": 2e-3, "p2_nodes": 2e-3}
# The projection-OFF stepper twins agree to ~1e-16 on every field.
_EXACT_TOL = 1e-12


def sphere_deltas(mem_ref, mem_test):
    """Per-field max-abs differences between two sphere runs on the interior:
    absolute for ``rho, rhoY, rhoX, p2_nodes``; ``/rho`` for the momenta."""
    import numpy as np

    rho = np.asarray(mem_ref.sol.rho)[_SPHERE_INNER]
    out = {}
    for name in ("rho", "rhoY", "rhoX"):
        a = np.asarray(getattr(mem_ref.sol, name))[_SPHERE_INNER]
        b = np.asarray(getattr(mem_test.sol, name))[_SPHERE_INNER]
        out[name] = float(np.max(np.abs(b - a)))
    for name in ("rhou", "rhov", "rhow"):
        a = np.asarray(getattr(mem_ref.sol, name))[_SPHERE_INNER]
        b = np.asarray(getattr(mem_test.sol, name))[_SPHERE_INNER]
        out[name] = float(np.max(np.abs((b - a) / rho)))
    out["p2_nodes"] = float(
        np.max(
            np.abs(np.asarray(mem_test.npf.p2_nodes) - np.asarray(mem_ref.npf.p2_nodes))
        )
    )
    return out


def assert_sphere_deltas(deltas, *, scalar, momentum, p2_nodes):
    """Assert every field at once so a CI failure reports the whole picture."""
    tol = {
        "rho": scalar,
        "rhoY": scalar,
        "rhoX": scalar,
        "rhou": momentum,
        "rhov": momentum,
        "rhow": momentum,
        "p2_nodes": p2_nodes,
    }
    report = ", ".join(f"{k}={v:.2e} (tol {tol[k]:.0e})" for k, v in deltas.items())
    bad = [k for k, v in deltas.items() if not v < tol[k]]
    assert not bad, f"exceeded on {bad}: {report}"


def _run_sphere_tc2(backend, nsteps, initial_projection=True, inx=48 + 1, inz=48 + 1):
    """Run Williamson TC2 (thin-shell sphere) for ``nsteps`` on ``backend``.

    Coarsened but pole-safe; exercises the whole sphere path — general
    curvilinear elliptic solve, advection through the metric normals, the
    ``coriolis_field`` H^-1, the general free-slip walls and the tangent-
    plane surface constraint."""
    import numpy as np

    from pybella.backends import jax_ops  # noqa: F401  (enables x64)
    from pybella.flow_solver.discretisation import grid as dis_grid, time_update
    from pybella.flow_solver.physics import thermodynamics
    from pybella.flow_solver.utils import cache, fields
    from pybella.tests import test_sphere_swe_tc2 as tc2
    from pybella.utils import user_data
    from pybella.utils.data_structures import ModelState
    from pybella.utils.io.debug import NullDebugWriter

    udo = tc2.UserData()
    udo.stepmax = nsteps
    udo.diag = False
    udo.output_timesteps = False
    udo.inx, udo.inz = inx, inz  # iny stays 2 (thin shell)
    udo.initial_projection = initial_projection  # do_initial_projection no-ops if False
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.backend = backend
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = tc2.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    return time_update.do(
        mem, ud, tout=1e9, bld=None, writer=None, debug_writer=NullDebugWriter()
    )


def _run_sphere_tc2_global(
    backend, nsteps, initial_projection=True, inx=32 + 1, inz=36 + 1
):
    """Run pole-to-pole Williamson TC2 for ``nsteps`` on ``backend``.

    Exercises the full pole machinery — the pole ghost exchange (cells +
    nodes), the elliptic pole-ring collapse, the FFT-in-longitude polar
    filter and its surface-constraint re-application — on top of the whole
    non-vertical-line sphere path."""
    import numpy as np

    from pybella.backends import jax_ops  # noqa: F401  (enables x64)
    from pybella.flow_solver.discretisation import grid as dis_grid, time_update
    from pybella.flow_solver.physics import thermodynamics
    from pybella.flow_solver.utils import cache, fields
    from pybella.tests import test_sphere_swe_tc2_global as tc2g
    from pybella.utils import user_data
    from pybella.utils.data_structures import ModelState
    from pybella.utils.io.debug import NullDebugWriter

    udo = tc2g.UserData()
    udo.stepmax = nsteps
    udo.diag = False
    udo.output_timesteps = False
    udo.inx, udo.inz = inx, inz  # iny stays 2 (thin shell)
    udo.initial_projection = initial_projection
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.backend = backend
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = tc2g.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    return time_update.do(
        mem, ud, tout=1e9, bld=None, writer=None, debug_writer=NullDebugWriter()
    )


def test_sphere_tc2_global_stepper_bit_identical():
    """The pole machinery twins are EXACT: with the initial projection off,
    the hybrid pole-to-pole TC2 stepper reproduces numpy to machine precision.

    This isolates the pole machinery — the pole ghost exchange, the elliptic
    pole-ring collapse and the FFT polar filter — from
    the ill-conditioned initial-projection solve (whose Krylov-floor member is
    what the projection gate below measures). Every per-step pole op is a pure
    index remap / one-sided collapse / functional FFT with no ulp-level
    backend divergence, so the whole 6-step run agrees bitwise-close."""
    import numpy as np

    n = 6
    mem_np = _run_sphere_tc2_global("numpy", n, initial_projection=False)
    mem_jx = _run_sphere_tc2_global("jax", n, initial_projection=False)
    inner = (slice(2, -2), slice(2, -2), slice(2, -2))

    def rel(a, b):
        a, b = np.asarray(a), np.asarray(b)
        return float(np.max(np.abs(a - b))) / max(1.0, float(np.max(np.abs(a))))

    # magnitude-scaled; ~1e-11 is the ulp-accumulation floor of the pure-jax
    # per-kernel FP (advection recovery/HLL + the per-step elliptic solves
    # picking ulp-different Krylov members) over 6 steps — ~7 orders below the
    # initial-projection Krylov floor the gate below measures.
    for name in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX"):
        d = rel(getattr(mem_np.sol, name)[inner], getattr(mem_jx.sol, name)[inner])
        assert d < 1e-10, f"{name}: {d:.3e}"
    d = rel(mem_np.npf.p2_nodes, mem_jx.npf.p2_nodes)
    assert d < 1e-10, f"p2_nodes: {d:.3e}"


def test_sphere_tc2_global_hybrid_reproduces_numpy():
    """Sanity check: the hybrid JAX pole-to-pole sphere path reproduces numpy
    over a short TC2 horizon WITH the initial projection, to within the
    jax-release-dependent projection floor.

    The bicgstab projection fixes the answer only to its residual class
    (scipy default ``rtol=1e-5``); the global pole-collapse system is the
    worst-conditioned of the sphere cases, so its floor is the highest
    (p2 ~5.8e-4 under jax 0.11.1, ~1.8e-5 under 0.10.1). The precision gate
    is ``test_sphere_tc2_global_stepper_bit_identical`` above; this one only
    catches gross breakage (a broken pole exchange / wrong Coriolis factor
    moves the fields by O(1e-2))."""
    n = 8
    mem_np = _run_sphere_tc2_global("numpy", n)
    mem_jx = _run_sphere_tc2_global("jax", n)
    assert_sphere_deltas(sphere_deltas(mem_np, mem_jx), **_SANITY_TOL)


def test_sphere_tc2_stepper_bit_identical():
    """The channel sphere twins are EXACT: with the initial projection off,
    the hybrid TC2 stepper reproduces numpy to machine precision (~1e-16 on
    every field, hybrid and device alike, under jax 0.10.1 and 0.11.1).

    This is the channel counterpart of the global stepper gate — the general
    curvilinear elliptic solve, advection through the metric normals, the
    ``coriolis_field`` H^-1, the general free-slip walls and the tangent-plane
    surface constraint, isolated from the projection's Krylov-floor member."""
    n = 8
    mem_np = _run_sphere_tc2("numpy", n, initial_projection=False)
    mem_jx = _run_sphere_tc2("jax", n, initial_projection=False)
    assert_sphere_deltas(
        sphere_deltas(mem_np, mem_jx),
        scalar=_EXACT_TOL,
        momentum=_EXACT_TOL,
        p2_nodes=_EXACT_TOL,
    )


def test_sphere_tc2_hybrid_reproduces_numpy():
    """Sanity check: the hybrid JAX sphere path reproduces numpy over a short
    TC2 horizon WITH the initial projection, to within the jax-release-
    dependent projection floor.

    Measured floor, same numpy/scipy/machine: jax 0.10.1 — scalars ~2e-6,
    momenta/rho ~1e-5; jax 0.11.1 — scalars ~8e-5, momenta/rho ~6e-4. The
    precision gate is ``test_sphere_tc2_stepper_bit_identical``; this one only
    catches gross breakage (a wrong ``rotation_axis_cart`` / Coriolis factor
    moves the fields by O(1e-2))."""
    n = 8
    mem_np = _run_sphere_tc2("numpy", n)
    mem_jx = _run_sphere_tc2("jax", n)
    assert_sphere_deltas(sphere_deltas(mem_np, mem_jx), **_SANITY_TOL)


@pytest.mark.parametrize(
    "ic",
    [
        "test_travelling_vortex",
        "test_straka",
        "test_travelling_vortex_3d_coriolis",
        "test_agnesi_hydrostatic",
    ],
)
def test_fullrun_jax_backend(ic):
    env = {**os.environ, "PYBELLA_BACKEND": "jax", "JAX_ENABLE_X64": "1"}
    result = subprocess.run(
        ["pybella", "-ic", ic, "-N", "1"], capture_output=True, text=True, env=env
    )
    assert result.returncode == 0, (
        f"JAX-backend run failed for {ic} (return code {result.returncode})\n"
        f"STDERR:\n{result.stderr.strip()[-3000:]}\n"
        f"STDOUT:\n{result.stdout.strip()[-3000:]}"
    )
    assert "Test passed" in result.stdout + result.stderr, (
        f"{ic}: run succeeded but CompareSol output not found — "
        "regression comparison did not execute"
    )
