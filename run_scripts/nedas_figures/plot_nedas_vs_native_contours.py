"""p2_nodes member CONTOUR panels at t=3.0: pyBELLA native LETKF vs NEDAS ETKF.
Generates the momentum-only and the Fig-10 all-quantities comparison."""

import glob
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MAIN = "/home/ray/git-projects/pybella/outputs/test_travelling_vortex"
WT = "/home/ray/git-projects/pybella-nedas/outputs"
T, K, INNER = 3.0, 10, (slice(2, -2), slice(2, -2))
NLEV = 12  # even count: no zero contour


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
    assert len(dirs) == K, (work, tag, len(dirs))
    return [np.load(d + "/p2_nodes.npy") for d in dirs]


truth = native_members(
    f"{WT}/test_travelling_vortex/test_travelling_vortex_ensemble=1_64_64_3.000000_truth_ib-0.h5",
    nmem=1,
)
noda = nedas_members("nedas_tv_noda", "prior")


def make(rows, title, out):
    fig, axes = plt.subplots(len(rows), K, figsize=(2.0 * K, 2.0 * len(rows)))
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
                    linewidths=0.7,
                    negative_linestyles="dashed",
                )
                ax.set_aspect("equal")
            else:
                ax.set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
            if c == 0:
                ax.set_ylabel("%s\n±%.2e" % (label, vmax), fontsize=9)
            if r == len(rows) - 1:
                ax.set_xlabel("mem %i" % c, fontsize=9)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


make(
    [
        ("truth", truth),
        ("noda\n(NEDAS fcst)", noda),
        (
            "EnDA\npyBELLA LETKF",
            native_members(
                f"{MAIN}/test_travelling_vortex_ensemble=10_64_64_3.000000_wdawloc_rhou_rhov_wda_ib-0.h5"
            ),
        ),
        ("EnDA\nNEDAS ETKF", nedas_members("nedas_tv_enda", "post")),
        (
            "EnDAB\npyBELLA LETKF",
            native_members(
                f"{MAIN}/test_travelling_vortex_ensemble=10_64_64_3.000000_wdawloc_rhou_rhov_wda_ib-0_cont_blend_fs=1_ts=0.h5"
            ),
        ),
        ("EnDAB\nNEDAS ETKF", nedas_members("nedas_tv_endab", "post")),
    ],
    "TV p2_nodes at t = 3.000, MOMENTUM-ONLY obs {rhou,rhov}, %i levels/row (dashed = negative) — native LETKF vs NEDAS ETKF"
    % NLEV,
    f"{WT}/nedas_vs_native_p2_t3_contours.png",
)
make(
    [
        ("truth", truth),
        ("noda\n(NEDAS fcst)", noda),
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
    ],
    "TV p2_nodes at t = 3.000, ALL-quantities obs (Fig-10 config), %i levels/row (dashed = negative) — native LETKF vs NEDAS ETKF"
    % NLEV,
    f"{WT}/nedas_vs_native_p2_t3_contours_all.png",
)
