"""Full regression runs on the device-resident JAX backend (phase B gate).

Complete production cases under ``PYBELLA_BACKEND=jax-device`` must pass
``CompareSol.test_do`` against the *same stored golden-master targets* as
the numpy and hybrid backends. The four CI cases cover 2D periodic
advection+elliptic, diffusion + x-WALLs, 3D elliptic + full Coriolis, and
terrain-following coordinates; the remaining cases were validated when the
phase landed — run them ad hoc with
``PYBELLA_BACKEND=jax-device pytest test_scripts/test_flow_solver.py``.

Skips cleanly when jax is not installed.
"""

import importlib.util
import os
import subprocess

import pytest

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("jax") is None, reason="jax not installed"
)


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
