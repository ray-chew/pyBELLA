"""Vmap gate: ensemble-vmapped vs looped batch forecasts through NEDAS.

Runs the travelling-vortex EnDA OSSE (2 cycles, K=10, initial blending off —
the jax-device guard rejects blending) twice on the SAME device with
`ens_run_strategy: batch`: once with `ens_batch_mode: vmap` (one vmapped
lockstep window integration over the member axis) and once with `loop` (the
identical batch-min-dt driver stepping members one by one). Both share the
dt sequence by construction, so any difference is pure vmap/jit reordering;
the gate is max-abs <= 1e-8 over every prior/post member checkpoint
(measured: bitwise identical on CPU jax and H100, jax 0.10.1).

Needs the regenerable TV obs file (run_scripts/osse_mwr2022.py tv --runs
obs); exits 0 with a SKIP message when it is missing. Backend: forces
PYBELLA_BACKEND=jax-device + JAX_ENABLE_X64=1 in the subprocesses (any
device jax finds — GPU on the cluster, CPU otherwise).

Run: python test_scripts/test_nedas_vmap_gate.py
"""

import glob
import os
import subprocess
import sys

import numpy as np
import yaml

OBS_GLOB = "outputs/test_travelling_vortex/*ensemble=1_64_64*_obs.h5"
BASE_CONFIG = "run_scripts/nedas_tv_osse.yml"
GATE = 1e-8
WORK = {mode: f"outputs/nedas_vmap_gate_{mode}" for mode in ("vmap", "loop")}


def build_config(mode: str) -> str:
    with open(BASE_CONFIG) as f:
        cfg = yaml.safe_load(f)
    cfg["work_dir"] = WORK[mode]
    cfg["time_end"] = "2001-01-01T00:30:00Z"  # 2 cycles: enough for a
    # forecast->analysis->forecast round trip through the batch driver
    model = cfg["model_def"]["pybella"]
    model["ens_run_strategy"] = "batch"
    model["ens_batch_mode"] = mode
    model["ud_overrides"]["initial_blending"] = False
    path = f"{WORK[mode]}_config.yml"
    with open(path, "w") as f:
        yaml.safe_dump(cfg, f)
    return path


def main() -> int:
    if not glob.glob(OBS_GLOB):
        print("SKIP: no TV obs file (run run_scripts/osse_mwr2022.py tv --runs obs)")
        return 0

    env = {
        **os.environ,
        "PYBELLA_BACKEND": "jax-device",
        "JAX_ENABLE_X64": "1",
    }
    for mode in ("vmap", "loop"):
        subprocess.run(["rm", "-rf", WORK[mode]], check=True)
        subprocess.run(
            [sys.executable, "run_scripts/nedas_run.py", "-c", build_config(mode)],
            check=True,
            env=env,
        )

    base = {m: os.path.join(WORK[m], "memory", "model", "pybella") for m in WORK}
    files = sorted(glob.glob(os.path.join(base["vmap"], "*", "*", "*.npy")))
    assert files, "no npy checkpoints written by the vmap run"
    worst, worst_at, checked = 0.0, "", 0
    for fv in files:
        rel = os.path.relpath(fv, base["vmap"])
        fl = os.path.join(base["loop"], rel)
        assert os.path.exists(fl), f"loop run missing checkpoint {rel}"
        d = float(np.max(np.abs(np.load(fv) - np.load(fl))))
        if d > worst:
            worst, worst_at = d, rel
        checked += 1
    assert worst <= GATE, (
        f"vmap vs loop max-abs {worst:.3e} at {worst_at} exceeds gate {GATE:.0e}"
    )
    print(
        f"PASS: {checked} member checkpoints, vmap vs loop max-abs "
        f"{worst:.3e} (gate {GATE:.0e})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
