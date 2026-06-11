"""Full regression runs on the JAX backend (components 2+3 end-to-end).

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
validated once on this backend when the port landed; run them ad hoc with
``PYBELLA_BACKEND=jax pytest test_scripts/test_flow_solver.py``.

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
