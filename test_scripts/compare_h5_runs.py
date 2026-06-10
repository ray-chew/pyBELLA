"""Bit-for-bit comparison of two pyBELLA output H5 files (Phase gate tool).

CompareSol's per-field tolerances (1e-5) are necessary but not sufficient for
the axial-agnosticity refactor's pure phases, which must be *bit-identical*
for the default vertical axis. This walks every dataset common to two run
files and reports max |a - b|; exit code 1 if any dataset differs (or is
missing from one side).

Usage:
    python test_scripts/compare_h5_runs.py baseline.h5 candidate.h5 [--tol 0]

Typical workflow: copy ./outputs/test_<case>/<case>_<N>_<M>.h5 for all cases
into a scratch baseline dir before starting a phase, rerun the suite after,
then compare pairwise.
"""

import argparse
import sys

import h5py
import numpy as np


def collect(h5):
    out = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            out[name] = obj[...]

    h5.visititems(visit)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline")
    ap.add_argument("candidate")
    ap.add_argument("--tol", type=float, default=0.0)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    with h5py.File(args.baseline, "r") as fa, h5py.File(args.candidate, "r") as fb:
        a, b = collect(fa), collect(fb)

    only_a = sorted(set(a) - set(b))
    only_b = sorted(set(b) - set(a))
    failures = []
    worst = []

    for name in sorted(set(a) & set(b)):
        if a[name].shape != b[name].shape:
            failures.append(f"{name}: shape {a[name].shape} vs {b[name].shape}")
            continue
        if a[name].dtype.kind not in "fiu":
            continue
        diff = np.max(np.abs(np.asarray(a[name], float) - np.asarray(b[name], float)))
        worst.append((diff, name))
        if diff > args.tol:
            failures.append(f"{name}: max|diff| = {diff:.3e}")

    if only_a:
        failures.append(
            f"only in baseline: {only_a[:5]}{'...' if len(only_a) > 5 else ''}"
        )
    if only_b:
        failures.append(
            f"only in candidate: {only_b[:5]}{'...' if len(only_b) > 5 else ''}"
        )

    if not args.quiet:
        worst.sort(reverse=True)
        print(f"{len(a)} / {len(b)} datasets; 5 largest diffs:")
        for d, n in worst[:5]:
            print(f"  {d:.3e}  {n}")

    if failures:
        print(f"FAIL ({len(failures)} issues, tol={args.tol:g}):")
        for f in failures[:20]:
            print(f"  {f}")
        sys.exit(1)
    print(f"OK: bit-identical at tol={args.tol:g}")


if __name__ == "__main__":
    main()
