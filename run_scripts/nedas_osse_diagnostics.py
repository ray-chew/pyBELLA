"""RMSE / spread diagnostics for NEDAS-driven pyBELLA OSSEs.

Reads the npy memory checkpoints written by a nedas_run.py experiment
(``<work_dir>/memory/model/pybella/<tstr>/<prior|post>_mem<NNN>/<attr>.npy``,
ghost-free (y, x) inner-domain arrays) and compares against the NATIVE truth
run HDF5, with the same conventions as run_scripts/osse_diagnostics.py:
ensemble-mean RMSE vs truth and grid-mean ensemble std (ddof=1), inner domain.
'prior' maps to the native '_before_da' state, 'post' to '_after_full_step'.

Usage:

    python run_scripts/nedas_osse_diagnostics.py \
        --work outputs/nedas_tv_osse \
        --truth outputs/test_travelling_vortex/..._truth_....h5 \
        --time-start 2001-01-01T00:00:00 --out outputs/nedas_diag_tv
"""

import argparse
import csv
import glob
import os
from datetime import datetime

import h5py
import numpy as np

def _inner(ndim):
    return (slice(2, -2),) * ndim  # gate-checked cases: 2 ghosts per axis


def read_truth(h5, field, t):
    label = "%s_ensemble_mem=0_%.3f_after_full_step" % (field, t)
    if label not in h5[field]:
        return None  # e.g. t=0: the truth run has no output before first tout
    data = np.squeeze(h5[field][label][:])
    return data[_inner(data.ndim)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True, help="NEDAS work_dir")
    ap.add_argument("--truth", required=True, help="native truth-run HDF5")
    ap.add_argument("--time-start", default="2001-01-01T00:00:00")
    ap.add_argument("--out", required=True, help="output directory")
    args = ap.parse_args()

    t0 = datetime.fromisoformat(args.time_start)
    base = os.path.join(args.work, "memory", "model", "pybella")
    os.makedirs(args.out, exist_ok=True)

    rows = []
    w_rows = []
    with h5py.File(args.truth, "r") as truth_h5:
        for tstr in sorted(os.listdir(base)):
            t = round(
                (datetime.strptime(tstr, "%Y%m%d_%H%M") - t0).total_seconds() / 3600.0,
                3,
            )
            for tag, state in (("prior", "before_da"), ("post", "after_full_step")):
                mem_dirs = sorted(glob.glob(os.path.join(base, tstr, f"{tag}_mem*")))
                if not mem_dirs:
                    continue
                fields = [
                    os.path.splitext(os.path.basename(p))[0]
                    for p in glob.glob(os.path.join(mem_dirs[0], "*.npy"))
                ]
                stacks = {}
                for field in sorted(fields):
                    members = np.array(
                        [np.load(os.path.join(d, field + ".npy")) for d in mem_dirs]
                    )
                    if members.ndim == 3:
                        # 2D snapshots are (y, x); native truth inner is (x, y)
                        members = members.transpose(0, 2, 1)
                    # 3D snapshots are already native (x, y, z) — no transpose
                    truth = read_truth(truth_h5, field, t)
                    if truth is None:
                        continue
                    stacks[field] = (members, truth)
                    mean = members.mean(axis=0)
                    rmse = np.sqrt(((mean - truth) ** 2).mean())
                    spread = members.std(axis=0, ddof=1).mean()
                    rows.append(
                        {
                            "field": field,
                            "t": t,
                            "state": state,
                            "rmse": rmse,
                            "spread": spread,
                            "K": len(mem_dirs),
                        }
                    )
                # w-probe (the hydrostatic-blending diagnostic, Phase E
                # headline): vertical velocity w = rhov/rho — DA-injected
                # imbalance shows as an acoustic w burst that EnDAB blends
                # away and EnDA does not
                if "rho" in stacks and "rhov" in stacks:
                    (rhos, rho_t) = stacks["rho"]
                    (rhovs, rhov_t) = stacks["rhov"]
                    w = rhovs / rhos
                    w_truth = rhov_t / rho_t
                    w_rows.append(
                        {
                            "t": t,
                            "state": state,
                            "rmse_w": np.sqrt(
                                ((w.mean(axis=0) - w_truth) ** 2).mean()
                            ),
                            "ens_mean_maxw": np.abs(w).max(
                                axis=tuple(range(1, w.ndim))
                            ).mean(),
                            "truth_maxw": np.abs(w_truth).max(),
                            "K": len(mem_dirs),
                        }
                    )

    csv_path = os.path.join(args.out, "rmse_spread.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["field", "t", "state", "rmse", "spread", "K"]
        )
        writer.writeheader()
        writer.writerows(rows)
    print("wrote %s (%i rows)" % (csv_path, len(rows)))

    if w_rows:
        w_path = os.path.join(args.out, "w_probe.csv")
        with open(w_path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "t", "state", "rmse_w", "ens_mean_maxw", "truth_maxw", "K",
                ],
            )
            writer.writeheader()
            writer.writerows(w_rows)
        print("wrote %s (%i rows)" % (w_path, len(w_rows)))
    for r in rows:
        print(
            "%-10s t=%-5.2f %-16s rmse=%.4e spread=%.4e"
            % (r["field"], r["t"], r["state"], r["rmse"], r["spread"])
        )


if __name__ == "__main__":
    main()
