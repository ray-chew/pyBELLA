"""Reproducibility-gate harness for behaviour-preserving refactors.

A refactor that is meant to change *no numerics* must leave every solver output
bit-for-bit unchanged. This harness makes that claim provable rather than
asserted:

    # before touching code, on a clean tree:
    python test_scripts/repro_gate.py capture --set fast

    # ... apply the refactor ...

    # after:
    python test_scripts/repro_gate.py check --set fast      # exit 0 == inert

``capture`` runs each case (``pybella -ic <case> -N 1``), records the per-field
``Max Abs`` gate line, and snapshots the output ``.h5`` to a scratch baseline
dir *outside* the repo. ``check`` reruns and, per case, requires (a) exit 0 (the
inline ``CompareSol`` gate stayed green) and (b) the output ``.h5`` is
bit-identical to the baseline at ``--tol 0`` via :mod:`compare_h5_runs`.
Bit-identity is the strong oracle: it implies the gate verdict and every
reported max-abs are unchanged.

This is the IRON RULE in code: a refactor must never regenerate a target,
loosen a tolerance, or change which fields are compared. The gates are the
oracle, not the edit.

Case sets
---------
fast    cheap non-terrain physics, run after every commit:
        travelling_vortex (2D elliptic), straka (diffusion + x-walls),
        travelling_vortex_3d_coriolis (3D + Coriolis).
terrain agnesi_hydrostatic, schaer_ridge — slow (per-process numba JIT compile
        of the 3D terrain kernels); run at phase boundaries only.
full    every golden-master regression case (fast + terrain + the rest).
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_BASELINE = os.path.expanduser("~/.cache/pybella-refactor/baseline")

CASE_SETS = {
    "fast": [
        "test_travelling_vortex",
        "test_straka",
        "test_travelling_vortex_3d_coriolis",
    ],
    "terrain": [
        "test_agnesi_hydrostatic",
        "test_schaer_ridge",
    ],
    "full": [
        "test_travelling_vortex",
        "test_travelling_vortex_3d_coriolis",
        "test_internal_long_wave",
        "test_igw_baldauf_brdar",
        "test_lamb_wave",
        "test_unstable_lamb",
        "test_swe_vortex",
        "test_straka",
        "test_agnesi_hydrostatic",
        "test_schaer_ridge",
        "test_blending_warm_bubble",
    ],
}


def _output_h5(case):
    """The final (non-``_old``) output h5 for a case, or None."""
    cands = [
        p
        for p in glob.glob(os.path.join(REPO, "outputs", case, f"{case}_*.h5"))
        if not p.endswith("_old.h5")
    ]
    return sorted(cands)[0] if cands else None


def _run(case, env=None):
    """Run a case; return (returncode, max_abs_lines)."""
    proc = subprocess.run(
        ["pybella", "-ic", case, "-N", "1"],
        capture_output=True,
        text=True,
        cwd=REPO,
        env=env,
    )
    # the inline CompareSol gate logs to stderr; scan both streams.
    max_abs = [
        ln.strip()
        for ln in (proc.stdout + proc.stderr).splitlines()
        if "Max Abs" in ln and "Test " in ln
    ]
    return proc.returncode, max_abs


def cmd_capture(cases, baseline):
    os.makedirs(baseline, exist_ok=True)
    bad = []
    for case in cases:
        rc, max_abs = _run(case)
        h5 = _output_h5(case)
        if rc != 0 or h5 is None:
            bad.append(case)
            print(f"  [SKIP] {case}: exit={rc}, h5={'missing' if h5 is None else 'ok'}")
            continue
        shutil.copy(h5, os.path.join(baseline, f"{case}.h5"))
        with open(os.path.join(baseline, f"{case}.maxabs.txt"), "w") as fh:
            fh.write("\n".join(max_abs) + "\n")
        print(f"  [OK]   {case}: {len(max_abs)} fields, h5 snapshot saved")
    if bad:
        print(f"\nFAILED to capture a green baseline for: {bad}")
        return 1
    print(f"\nbaseline captured in {baseline}")
    return 0


def cmd_check(cases, baseline, tol):
    fails = []
    for case in cases:
        ref_h5 = os.path.join(baseline, f"{case}.h5")
        if not os.path.exists(ref_h5):
            print(f"  [MISS] {case}: no baseline; run `capture` first")
            fails.append(case)
            continue
        rc, _ = _run(case)
        if rc != 0:
            print(f"  [GATE] {case}: pybella exit={rc} (CompareSol gate FAILED)")
            fails.append(case)
            continue
        cand_h5 = _output_h5(case)
        diff = subprocess.run(
            ["python", os.path.join(HERE, "compare_h5_runs.py"), ref_h5, cand_h5,
             "--tol", str(tol), "--quiet"],
            capture_output=True,
            text=True,
            cwd=REPO,
        )
        if diff.returncode != 0:
            print(f"  [DIFF] {case}: NOT bit-identical:\n{diff.stdout.strip()}")
            fails.append(case)
        else:
            print(f"  [OK]   {case}: green + bit-identical at tol={tol:g}")
    if fails:
        print(f"\nREPRO GATE FAILED for: {fails}")
        return 1
    print(f"\nrepro gate PASSED ({len(cases)} cases inert)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=["capture", "check"])
    ap.add_argument("--set", dest="set_name", default="fast", choices=list(CASE_SETS))
    ap.add_argument("--cases", nargs="*", help="explicit case list (overrides --set)")
    ap.add_argument("--baseline", default=DEFAULT_BASELINE)
    ap.add_argument("--tol", type=float, default=0.0)
    args = ap.parse_args()

    cases = args.cases or CASE_SETS[args.set_name]
    if args.mode == "capture":
        sys.exit(cmd_capture(cases, args.baseline))
    sys.exit(cmd_check(cases, args.baseline, args.tol))


if __name__ == "__main__":
    main()
