"""Member-gallery quick-look: vertical (x,y) slices of w = rhov/rho at mid-z,
truth + all 10 members, EnDA vs EnDAB side by side, shared diverging scale.

Usage: python plot_members_w.py <t_nondim> <prior|post> <outfile.png>
"""

import os
import sys
from datetime import datetime, timedelta

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

T = float(sys.argv[1])
STATE = sys.argv[2] if len(sys.argv) > 2 else "prior"
OUT = sys.argv[3] if len(sys.argv) > 3 else f"outputs/e_figs/w_members_t{int(T)}_{STATE}.png"

TRUTH = sorted(
    __import__("glob").glob(
        "outputs/test_internal_long_wave/*ensemble=1_65_16_16*truth_ib-0.h5"
    )
)[-1]
K = 10
ZI = 8  # mid-z slice
XLIM, YLIM = (-300.0, 300.0), (0.0, 1.0)
W_SCALE = 1e6  # display in 1e-6 nondim units

t0 = datetime(2001, 1, 1)
tstr = (t0 + timedelta(hours=T)).strftime("%Y%m%d_%H%M")


def member_w(work, m):
    d = f"{work}/memory/model/pybella/{tstr}/{STATE}_mem{m:03d}"
    rho = np.load(f"{d}/rho.npy")
    rhov = np.load(f"{d}/rhov.npy")
    return (rhov / rho)[:, :, ZI].T * W_SCALE  # (y, x) for plotting


def truth_w():
    with h5py.File(TRUTH, "r") as f:
        lab = "%s_ensemble_mem=0_%.3f_after_full_step"
        rho = np.squeeze(f["rho"][lab % ("rho", T)][:])[2:-2, 2:-2, 2:-2]
        rhov = np.squeeze(f["rhov"][lab % ("rhov", T)][:])[2:-2, 2:-2, 2:-2]
    return (rhov / rho)[:, :, ZI].T * W_SCALE


runs = {"EnDA": "outputs/nedas_igw3d_enda", "EnDAB": "outputs/nedas_igw3d_endab"}
wt = truth_w()
fields = {
    lab: [wt] + [member_w(work, m) for m in range(1, K + 1)]
    for lab, work in runs.items()
}

# shared symmetric scale (robust: 99.5th pct over everything)
lim = max(
    np.percentile(np.abs(np.array(fields[lab])), 99.5) for lab in fields
)

nrows = K + 1
fig, axes = plt.subplots(
    nrows, 2, figsize=(11.5, 0.62 * nrows + 1.2), sharex=True, sharey=True,
    constrained_layout=True,
)
x = np.linspace(XLIM[0], XLIM[1], wt.shape[1] + 1)
y = np.linspace(YLIM[0], YLIM[1], wt.shape[0] + 1)
row_labels = ["truth"] + [f"m{m}" for m in range(1, K + 1)]

for j, lab in enumerate(runs):
    axes[0, j].set_title(f"{lab}  ({STATE}, t = {T:g})", fontsize=11)
    for i in range(nrows):
        ax = axes[i, j]
        pm = ax.pcolormesh(
            x, y, fields[lab][i], cmap="RdBu_r", vmin=-lim, vmax=lim,
            rasterized=True,
        )
        ax.tick_params(length=2, labelsize=7)
        if j == 0:
            ax.set_ylabel(row_labels[i], rotation=0, ha="right", va="center",
                          fontsize=9, labelpad=12)
for ax in axes[-1]:
    ax.set_xlabel("x (nondim)", fontsize=9)

cb = fig.colorbar(pm, ax=axes, shrink=0.5, pad=0.01, aspect=40)
cb.set_label(r"w = rhov/rho  [$10^{-6}$ nondim]", fontsize=9)
fig.suptitle(
    f"igw3d dev 65x16x16, K=10 — vertical (x,y) slices at mid-z (k={ZI}); "
    "top row = truth (identical in both columns)",
    fontsize=10,
)
if os.path.dirname(OUT):
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=150)
print("wrote", OUT, f"(color scale ±{lim:.2f}e-6)")
