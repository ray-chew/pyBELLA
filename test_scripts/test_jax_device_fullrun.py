"""Full regression runs on the device-resident JAX backend (phase B gate).

Complete production cases under ``PYBELLA_BACKEND=jax-device`` must pass
``CompareSol.test_do`` against the *same stored golden-master targets* as
the numpy and hybrid backends. The four CI cases cover 2D periodic
advection+elliptic, diffusion + x-WALLs, 3D elliptic + full Coriolis, and
terrain-following coordinates; the remaining cases were validated when the
phase landed — run them ad hoc with
``PYBELLA_BACKEND=jax-device pytest test_scripts/test_flow_solver.py``.

The sphere path (non-vertical-line metric: general e_up buoyancy, general
H^-1, ``coriolis_field``, the well-balanced gravity fill and free-slip
walls, the tangent-plane surface constraint) is gated in-process by
``test_sphere_tc2_device_reproduces_numpy`` — a short Williamson TC2 run
compared jax-device-vs-numpy, not against the stored target, because the
initial-projection Krylov floor (~2e-5 in the momenta) sits above the 1e-5
regression tolerance (see dev_notes/sphere.md, HARD-WON pt 2).

Skips cleanly when jax is not installed.
"""

import importlib.util
import os
import subprocess

import pytest

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("jax") is None, reason="jax not installed"
)


def test_sphere_tc2_device_reproduces_numpy():
    """The device-resident JAX sphere path reproduces numpy over a short TC2
    horizon, at the documented initial-projection Krylov floor (same
    argument as the hybrid gate in test_jax_fullrun.py). Exercises the
    on-device general e_up buoyancy, general H^-1, coriolis_field, general
    walls and the surface constraint end to end."""
    import numpy as np

    from test_jax_fullrun import _run_sphere_tc2

    n = 8
    mem_np = _run_sphere_tc2("numpy", n)
    mem_dv = _run_sphere_tc2("jax-device", n)
    assert getattr(mem_dv, "_device_compile_count", None) == 2  # parity 0/1
    inner = (slice(2, -2), slice(2, -2), slice(2, -2))
    rho = mem_np.sol.rho[inner]

    for name in ("rho", "rhoY", "rhoX"):
        d = float(
            np.max(
                np.abs(
                    getattr(mem_dv.sol, name)[inner] - getattr(mem_np.sol, name)[inner]
                )
            )
        )
        assert d < 1e-5, f"{name} (absolute): {d:.3e}"
    for name in ("rhou", "rhov", "rhow"):
        d = float(
            np.max(
                np.abs(
                    (
                        getattr(mem_dv.sol, name)[inner]
                        - getattr(mem_np.sol, name)[inner]
                    )
                    / rho
                )
            )
        )
        assert d < 6e-5, f"{name} (momentum/rho): {d:.3e}"
    d = float(np.max(np.abs(mem_dv.npf.p2_nodes - mem_np.npf.p2_nodes)))
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
