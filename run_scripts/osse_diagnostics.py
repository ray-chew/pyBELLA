"""RMSE / spread diagnostics for the MWR-2022 OSSE reproductions.

For each requested time and field, computes on the inner domain:

* ensemble-mean RMSE vs the truth run,
* ensemble spread (grid-mean ensemble standard deviation),

for both the forecast (``_before_da``, where present) and the analysis
(``_after_full_step``) states, writes a CSV of the numbers and a PNG of the
RMSE/spread time series per field.

Usage (from the repository root):

    python run_scripts/osse_diagnostics.py \
        --truth outputs/.../..._ensemble=1_..._truth_....h5 \
        --ens   outputs/.../..._ensemble=10_..._wda....h5 \
        --times 0.25 0.5 ... --members 10 --out outputs/osse_diag_tv_wda
"""

import argparse
import csv
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np

FIELDS = ["rho", "rhou", "rhov", "rhoY", "p2_nodes"]
# ghost widths of the gate-checked cases (igx = igy = 2)
INNER = (slice(2, -2), slice(2, -2))


def read_field(h5, field, member, t, suffix):
    label = "%s_ensemble_mem=%i_%.3f_%s" % (field, member, t, suffix)
    data = h5[field][label][:]
    return np.squeeze(data)[INNER]


def rmse_and_spread(ens_h5, truth_h5, field, N, t, suffix):
    truth = read_field(truth_h5, field, 0, t, "after_full_step")
    members = np.array([read_field(ens_h5, field, n, t, suffix) for n in range(N)])
    mean = members.mean(axis=0)
    rmse = np.sqrt(((mean - truth) ** 2).mean())
    spread = members.std(axis=0, ddof=1).mean()
    return rmse, spread


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth", required=True, help="truth-run HDF5 file")
    ap.add_argument("--ens", required=True, help="ensemble-run HDF5 file")
    ap.add_argument("--times", required=True, nargs="+", type=float)
    ap.add_argument("--members", type=int, default=10)
    ap.add_argument("--fields", nargs="*", default=FIELDS)
    ap.add_argument("--out", required=True, help="output directory")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rows = []
    with h5py.File(args.truth, "r") as truth_h5, h5py.File(args.ens, "r") as ens_h5:
        for field in args.fields:
            for t in args.times:
                for suffix in ("before_da", "after_full_step"):
                    label = "%s_ensemble_mem=0_%.3f_%s" % (field, t, suffix)
                    if label not in ens_h5[field]:
                        continue  # e.g. no _before_da in runs without DA
                    rmse, spread = rmse_and_spread(
                        ens_h5, truth_h5, field, args.members, t, suffix
                    )
                    rows.append(
                        {
                            "field": field,
                            "t": t,
                            "state": suffix,
                            "rmse": rmse,
                            "spread": spread,
                        }
                    )

    csv_path = os.path.join(args.out, "rmse_spread.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["field", "t", "state", "rmse", "spread"]
        )
        writer.writeheader()
        writer.writerows(rows)
    print("wrote %s (%i rows)" % (csv_path, len(rows)))

    for field in args.fields:
        frows = [r for r in rows if r["field"] == field]
        if not frows:
            continue
        fig, ax = plt.subplots(figsize=(7, 4))
        for suffix, marker in (("before_da", "x"), ("after_full_step", "o")):
            srows = [r for r in frows if r["state"] == suffix]
            if srows:
                ax.plot(
                    [r["t"] for r in srows],
                    [r["rmse"] for r in srows],
                    marker=marker,
                    label="RMSE (%s)" % suffix,
                )
                ax.plot(
                    [r["t"] for r in srows],
                    [r["spread"] for r in srows],
                    marker=marker,
                    linestyle="--",
                    label="spread (%s)" % suffix,
                )
        ax.set_xlabel("t")
        ax.set_yscale("log")
        ax.set_title(field)
        ax.legend(fontsize=8)
        png_path = os.path.join(args.out, "rmse_spread_%s.png" % field)
        fig.savefig(png_path, bbox_inches="tight", dpi=150)
        plt.close(fig)
        print("wrote %s" % png_path)


if __name__ == "__main__":
    main()
