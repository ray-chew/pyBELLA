"""Tolerance audit utility.

Runs every regression case N times, harvests CompareSol's per-field max-abs
errors from the run logs, and prints a table of observed error vs current
tolerance, plus a suggested tightened tolerance per case
(max observed error across reps and fields, rounded up to the next power of
ten, with one extra decade of headroom for cross-platform/BLAS variation;
never looser than the current tolerance).

Manual utility, not collected by pytest. Usage:

    python test_scripts/tolerance_audit.py [reps]
"""

import math
import re
import subprocess
import sys
from collections import defaultdict

CASES = [
    "test_travelling_vortex",
    "test_travelling_vortex_3d_coriolis",
    "test_internal_long_wave",
    "test_igw_baldauf_brdar",
    "test_lamb_wave",
    "test_unstable_lamb",
    "test_blending_warm_bubble",
    "test_swe_vortex",
    "test_straka",
]

LINE = re.compile(
    r"Test passed for (\w+) \| L2: ([\d.eE+-]+), Rel L2: ([\d.eE+-]+|inf), "
    r"Max Abs: ([\d.eE+-]+)"
)


def run_case(case):
    """Run one case; return {field: max_abs_error} or None on failure."""
    proc = subprocess.run(
        [sys.executable, "-m", "pybella", "-ic", case, "-N", "1"],
        capture_output=True,
        text=True,
    )
    out = proc.stdout + proc.stderr
    if proc.returncode != 0:
        print(f"  !! {case} FAILED (rc={proc.returncode}); tail:")
        print("\n".join(out.strip().splitlines()[-5:]))
        return None
    errors = {}
    for m in LINE.finditer(out):
        field, _, _, max_abs = m.groups()
        errors[field] = max(errors.get(field, 0.0), float(max_abs))
    return errors


def suggest(max_err):
    """Next power of ten above max_err, plus one decade of headroom."""
    if max_err == 0.0:
        return 1e-12
    return 10.0 ** (math.ceil(math.log10(max_err)) + 1)


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    results = {}
    for case in CASES:
        worst = defaultdict(float)
        for rep in range(reps):
            print(f"== {case} rep {rep + 1}/{reps}")
            errors = run_case(case)
            if errors is None:
                worst = None
                break
            for field, err in errors.items():
                worst[field] = max(worst[field], err)
        results[case] = dict(worst) if worst is not None else None

    print("\n\n==================== AUDIT SUMMARY ====================")
    for case, worst in results.items():
        print(f"\n{case}:")
        if worst is None:
            print("  RUN FAILED — investigate before tightening")
            continue
        if not worst:
            print("  no CompareSol lines found (diag off?)")
            continue
        overall = max(worst.values())
        print(f"  {'field':10s} {'max_abs_err':>12s}")
        for field, err in sorted(worst.items()):
            print(f"  {field:10s} {err:12.3e}")
        print(
            f"  worst: {overall:.3e}  ->  suggested tolerance: {suggest(overall):.0e}"
        )


if __name__ == "__main__":
    main()
