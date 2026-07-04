"""RB p2_nodes member CONTOUR panels at t=1.0: native LETKF vs NEDAS ETKF."""

import glob
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MAIN = "/home/ray/git-projects/pybella/outputs/output_rising_bubble"
WT = "/home/ray/git-projects/pybella-nedas/outputs"
T, K, INNER = 1.0, 10, (slice(2, -2), slice(2, -2))
NLEV = 12


def native_members(h5path, suffix="after_full_step", nmem=K):
    with h5py.File(h5path, "r") as f:
        return [
            np.squeeze(
                f["p2_nodes"]["p2_nodes_ensemble_mem=%i_%.3f_%s" % (n, T, suffix)][:]
            )[INNER].T
            for n in range(nmem)
        ]


def nedas_members(work, tag):
    dirs = sorted(
        glob.glob(f"{WT}/{work}/memory/model/pybella/20010101_0100/{tag}_mem*")
    )
    assert len(dirs) == K, (work, tag, len(dirs))
    return [np.load(d + "/p2_nodes.npy") for d in dirs]


rows = [
    (
        "truth",
        native_members(
            f"{WT}/output_rising_bubble/output_rising_bubble_ensemble=1_160_80_1.000000_CFLfixed_truth_ib-0.h5",
            nmem=1,
        ),
    ),
    (
        "noda\n(native)",
        native_members(
            f"{MAIN}/output_rising_bubble_ensemble=10_160_80_1.000000_CFLfixed_noda_ib-0.h5"
        ),
    ),
    (
        "EnDA\npyBELLA LETKF",
        native_members(
            f"{MAIN}/output_rising_bubble_ensemble=10_160_80_1.000000_wdawloc_rhou_rhov_CFLfixed_wda_ib-0.h5"
        ),
    ),
    ("EnDA\nNEDAS ETKF", nedas_members("nedas_bubble_enda", "post")),
    (
        "EnDAB\npyBELLA LETKF",
        native_members(
            f"{MAIN}/output_rising_bubble_ensemble=10_160_80_1.000000_wdawloc_rhou_rhov_CFLfixed_wda_ib-0_cont_blend_fs=1_ts=0.h5"
        ),
    ),
    ("EnDAB\nNEDAS ETKF", nedas_members("nedas_bubble_endab", "post")),
]

fig, axes = plt.subplots(len(rows), K, figsize=(2.4 * K, 1.5 * len(rows)))
for r, (label, mems) in enumerate(rows):
    vmax = max(np.abs(m).max() for m in mems)
    levels = np.linspace(-vmax, vmax, NLEV)
    for c in range(K):
        ax = axes[r, c]
        if c < len(mems):
            f = mems[c]
            y, x = np.mgrid[0 : f.shape[0], 0 : f.shape[1]]
            ax.contour(
                x,
                y,
                f,
                levels=levels,
                cmap="RdBu_r",
                vmin=-vmax,
                vmax=vmax,
                linewidths=0.6,
                negative_linestyles="dashed",
            )
            ax.set_aspect("equal")
        else:
            ax.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        if c == 0:
            ax.set_ylabel("%s\n±%.2e" % (label, vmax), fontsize=8)
        if r == len(rows) - 1:
            ax.set_xlabel("mem %i" % c, fontsize=8)
fig.suptitle(
    "RB p2_nodes at t = 1.000 (analysis state), momentum obs, %i levels/row (dashed = negative) — native LETKF vs NEDAS ETKF"
    % NLEV,
    fontsize=12,
)
fig.tight_layout(rect=[0, 0, 1, 0.965])
out = f"{WT}/nedas_vs_native_bubble_p2_t1_contours.png"
fig.savefig(out, dpi=110, bbox_inches="tight")
print("wrote", out)
