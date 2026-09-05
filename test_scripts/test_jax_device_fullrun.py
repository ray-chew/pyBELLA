"""Full regression runs on the device-resident JAX backend.

Complete production cases under ``PYBELLA_BACKEND=jax-device`` must pass
``CompareSol.test_do`` against the *same stored golden-master targets* as
the numpy and hybrid backends. The four CI cases cover 2D periodic
advection+elliptic, diffusion + x-WALLs, 3D elliptic + full Coriolis, and
terrain-following coordinates; the remaining cases were validated once on
this backend — run them ad hoc with
``PYBELLA_BACKEND=jax-device pytest test_scripts/test_flow_solver.py``.

The sphere path (non-vertical-line metric: general e_up buoyancy, general
H^-1, ``coriolis_field``, the well-balanced gravity fill and free-slip
walls, the tangent-plane surface constraint) is gated in-process by short
Williamson TC2 runs compared jax-device-vs-numpy, not against the stored
target. As in test_jax_fullrun.py the ``*_stepper_bit_identical`` gates
(initial projection OFF) are exact (~1e-16) and jax-release independent,
while the projection-ON ``*_reproduces_numpy`` gate is a loose sanity check
at the jax-release-dependent projection floor (see that module's docstring).

Skips cleanly when jax is not installed.
"""

import importlib.util
import os
import subprocess

import pytest

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("jax") is None, reason="jax not installed"
)


def test_sphere_tc2_device_stepper_bit_identical():
    """The device-resident channel sphere twins are EXACT: with the initial
    projection off, the jax-device TC2 stepper reproduces numpy to machine
    precision (~1e-16 on every field). Exercises the on-device general e_up
    buoyancy, general H^-1, coriolis_field, general walls and the surface
    constraint end to end, isolated from the projection's Krylov-floor
    member. This is the precision gate for the device channel path."""
    from test_jax_fullrun import _EXACT_TOL, _run_sphere_tc2, assert_sphere_deltas
    from test_jax_fullrun import sphere_deltas

    n = 8
    mem_np = _run_sphere_tc2("numpy", n, initial_projection=False)
    mem_dv = _run_sphere_tc2("jax-device", n, initial_projection=False)
    assert getattr(mem_dv, "_device_compile_count", None) == 2  # parity 0/1
    assert_sphere_deltas(
        sphere_deltas(mem_np, mem_dv),
        scalar=_EXACT_TOL,
        momentum=_EXACT_TOL,
        p2_nodes=_EXACT_TOL,
    )


def test_sphere_tc2_device_reproduces_numpy():
    """Sanity check: the device-resident JAX sphere path reproduces numpy over
    a short TC2 horizon WITH the initial projection, to within the
    jax-release-dependent projection floor (identical to the hybrid floor —
    the device and hybrid paths land on the same Krylov member). The precision
    gate is ``test_sphere_tc2_device_stepper_bit_identical``."""
    from test_jax_fullrun import _SANITY_TOL, _run_sphere_tc2, assert_sphere_deltas
    from test_jax_fullrun import sphere_deltas

    n = 8
    mem_np = _run_sphere_tc2("numpy", n)
    mem_dv = _run_sphere_tc2("jax-device", n)
    assert getattr(mem_dv, "_device_compile_count", None) == 2  # parity 0/1
    assert_sphere_deltas(sphere_deltas(mem_np, mem_dv), **_SANITY_TOL)


def test_sphere_tc2_global_device_reproduces_numpy():
    """The device-resident JAX pole-to-pole sphere path reproduces numpy over
    a short TC2 horizon. Exercises the on-device pole ghost fold
    (cells + nodes), the elliptic pole-ring collapse, and the FFT
    polar filter + surface constraint inside the device step loop, on top of
    the whole non-vertical-line sphere path — with the initial projection off,
    so the only floor is the pure-jax per-kernel ulp accumulation (the
    projection Krylov floor is measured by the hybrid gate in
    test_jax_fullrun.py); the device reproduces numpy essentially bitwise."""
    import numpy as np

    from test_jax_fullrun import _run_sphere_tc2_global

    n = 6
    mem_np = _run_sphere_tc2_global("numpy", n, initial_projection=False)
    mem_dv = _run_sphere_tc2_global("jax-device", n, initial_projection=False)
    assert getattr(mem_dv, "_device_compile_count", None) == 2  # parity 0/1
    inner = (slice(2, -2), slice(2, -2), slice(2, -2))

    def rel(a, b):
        a, b = np.asarray(a), np.asarray(b)
        return float(np.max(np.abs(a - b))) / max(1.0, float(np.max(np.abs(a))))

    for name in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX"):
        d = rel(getattr(mem_np.sol, name)[inner], getattr(mem_dv.sol, name)[inner])
        assert d < 1e-9, f"{name}: {d:.3e}"
    d = rel(mem_np.npf.p2_nodes, mem_dv.npf.p2_nodes)
    assert d < 1e-9, f"p2_nodes: {d:.3e}"


@pytest.mark.parametrize(
    "ic",
    [
        "test_travelling_vortex",
        "test_straka",
        "test_travelling_vortex_3d_coriolis",
        "test_agnesi_hydrostatic",
    ],
)
def test_fullrun_jax_device_backend(ic):
    env = {**os.environ, "PYBELLA_BACKEND": "jax-device", "JAX_ENABLE_X64": "1"}
    result = subprocess.run(
        ["pybella", "-ic", ic, "-N", "1"], capture_output=True, text=True, env=env
    )
    assert result.returncode == 0, (
        f"jax-device run failed for {ic} (return code {result.returncode})\n"
        f"STDERR:\n{result.stderr.strip()[-3000:]}\n"
        f"STDOUT:\n{result.stdout.strip()[-3000:]}"
    )
    out = result.stdout + result.stderr
    assert "Test passed" in out, f"{ic}: run succeeded but CompareSol output not found"
    assert (
        "device step" in out
    ), f"{ic}: run succeeded but the device-resident loop was not used"
