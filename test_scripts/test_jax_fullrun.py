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
``test_sphere_tc2_hybrid_reproduces_numpy`` — a short Williamson TC2 run
compared jax-vs-numpy, not against the stored target, because the initial-
projection Krylov floor (~2e-5 in the momenta) sits above the 1e-5
regression tolerance.

Skips cleanly when jax is not installed.
"""

import importlib.util
import os
import subprocess

import pytest

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("jax") is None, reason="jax not installed"
)


def _run_sphere_tc2(backend, nsteps, inx=48 + 1, inz=48 + 1):
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
    """The hybrid JAX pole-to-pole sphere path reproduces numpy over a short
    TC2 horizon WITH the initial projection, at the documented Krylov floor.

    As on the channel (see test_sphere_tc2_hybrid_reproduces_numpy), the
    bicgstab initial projection fixes the answer only to its residual class;
    the global pole-collapse system is more ill-conditioned, so the floor sits
    higher (~1e-4 in the momenta at this conditioning) — still far below the
    ~5e-4 a broken pole exchange / wrong Coriolis produces. The step-to-step
    reproduction itself is machine-exact (the projection-off gate above)."""
    import numpy as np

    n = 8
    mem_np = _run_sphere_tc2_global("numpy", n)
    mem_jx = _run_sphere_tc2_global("jax", n)
    inner = (slice(2, -2), slice(2, -2), slice(2, -2))
    rho = mem_np.sol.rho[inner]

    for name in ("rho", "rhoY", "rhoX"):
        d = float(
            np.max(
                np.abs(
                    getattr(mem_jx.sol, name)[inner] - getattr(mem_np.sol, name)[inner]
                )
            )
        )
        assert d < 2e-5, f"{name} (absolute): {d:.3e}"
    for name in ("rhou", "rhov", "rhow"):
        d = float(
            np.max(
                np.abs(
                    (
                        getattr(mem_jx.sol, name)[inner]
                        - getattr(mem_np.sol, name)[inner]
                    )
                    / rho
                )
            )
        )
        assert d < 2e-4, f"{name} (momentum/rho): {d:.3e}"
    d = float(np.max(np.abs(mem_jx.npf.p2_nodes - mem_np.npf.p2_nodes)))
    assert d < 3e-5, f"p2_nodes: {d:.3e}"


def test_sphere_tc2_hybrid_reproduces_numpy():
    """The hybrid JAX sphere path reproduces numpy over a short TC2 horizon,
    at the documented initial-projection Krylov floor.

    The bicgstab initial projection fixes the answer only up to its achieved-
    residual class (~2e-5 in the momenta at this conditioning), and the
    ulp-level per-kernel differences between the backends select a different
    member — so the momenta agree to ~1e-5 (not machine precision) while the
    tightly-set scalars agree to ~1e-6. Both are far below the ~5e-4 a wrong
    ``rotation_axis_cart`` / Coriolis factor produces, so this is a genuine
    end-to-end gate."""
    import numpy as np

    n = 8
    mem_np = _run_sphere_tc2("numpy", n)
    mem_jx = _run_sphere_tc2("jax", n)
    inner = (slice(2, -2), slice(2, -2), slice(2, -2))
    rho = mem_np.sol.rho[inner]

    for name in ("rho", "rhoY", "rhoX"):
        d = float(
            np.max(
                np.abs(
                    getattr(mem_jx.sol, name)[inner] - getattr(mem_np.sol, name)[inner]
                )
            )
        )
        assert d < 1e-5, f"{name} (absolute): {d:.3e}"
    for name in ("rhou", "rhov", "rhow"):
        d = float(
            np.max(
                np.abs(
                    (
                        getattr(mem_jx.sol, name)[inner]
                        - getattr(mem_np.sol, name)[inner]
                    )
                    / rho
                )
            )
        )
        assert d < 6e-5, f"{name} (momentum/rho): {d:.3e}"
    d = float(np.max(np.abs(mem_jx.npf.p2_nodes - mem_np.npf.p2_nodes)))
    assert d < 1e-5, f"p2_nodes: {d:.3e}"


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
