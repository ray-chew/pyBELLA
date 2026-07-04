"""All-quantities EnDAB/EnDA comparison, imshow variant."""

import glob
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MAIN = "/home/ray/git-projects/pybella/outputs/test_travelling_vortex"
WT = "/home/ray/git-projects/pybella-nedas/outputs"
T, K, INNER = 3.0, 10, (slice(2, -2), slice(2, -2))


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
        glob.glob(f"{WT}/{work}/memory/model/pybella/20010101_0300/{tag}_mem*")
    )
    return [np.load(d + "/p2_nodes.npy") for d in dirs]


truth = native_members(
    f"{WT}/test_travelling_vortex/test_travelling_vortex_ensemble=1_64_64_3.000000_truth_ib-0.h5",
    nmem=1,
)
rows = [
    ("truth", truth),
    ("noda\n(NEDAS fcst)", nedas_members("nedas_tv_noda", "prior")),
    (
        "EnDA(all)\npyBELLA LETKF",
        native_members(
            f"{MAIN}/test_travelling_vortex_ensemble=10_64_64_3.000000_wdawloc_all_wda_all_ib-0.h5"
        ),
    ),
    ("EnDA(all)\nNEDAS ETKF", nedas_members("nedas_tv_enda_all", "post")),
    (
        "EnDAB(all)\npyBELLA LETKF",
        native_members(
            f"{MAIN}/test_travelling_vortex_ensemble=10_64_64_3.000000_wdawloc_all_wda_all_ib-0_cont_blend_fs=1_ts=0.h5"
        ),
    ),
    ("EnDAB(all)\nNEDAS ETKF", nedas_members("nedas_tv_endab_all", "post")),
]

fig, axes = plt.subplots(len(rows), K, figsize=(2.0 * K, 2.0 * len(rows)))
for r, (label, mems) in enumerate(rows):
    vmax = max(np.abs(m).max() for m in mems)
    for c in range(K):
        ax = axes[r, c]
        if c < len(mems):
            ax.imshow(mems[c], origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        else:
            ax.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        if c == 0:
            ax.set_ylabel("%s\n±%.2e" % (label, vmax), fontsize=9)
        if r == len(rows) - 1:
            ax.set_xlabel("mem %i" % c, fontsize=9)
fig.suptitle(
    "TV p2_nodes at t = 3.000, ALL-quantities obs (Fig-10 config), rows scaled per run — native LETKF vs NEDAS ETKF",
    fontsize=12,
)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = f"{WT}/nedas_vs_native_p2_t3_imshow_all.png"
fig.savefig(out, dpi=110, bbox_inches="tight")
print("wrote", out)
