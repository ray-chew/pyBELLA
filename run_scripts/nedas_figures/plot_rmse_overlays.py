"""RMSE/spread vs time overlays incl. t=0: native LETKF vs NEDAS ETKF.

t=0: the initial ensemble is IDENTICAL across engines and variants (same
seed-888 chain), and the truth IC regenerates from sol_init (truth aux) —
so every curve starts from one common initial-error point. For the bubble,
the NEDAS pre-DA spinup priors (t=0.05..0.45, forecast states) are added as
a dotted line to show what initial blending does before the first analysis.
"""

import csv, glob, importlib, os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

WT = "/home/ray/git-projects/pybella-nedas/outputs"
NAT = "/home/ray/git-projects/pybella/outputs"


def truth_ic(ic_key, aux):
    from pybella.utils import user_data
    from pybella.flow_solver.discretisation import grid as dis_grid
    from pybella.flow_solver.physics import thermodynamics as gd_th
    from pybella.flow_solver.utils import fields
    from pybella.interfaces.ic_config import IC_MODULES

    mod = importlib.import_module(IC_MODULES[ic_key])
    ud = user_data.UserDataInit(**vars(mod.UserData()))
    ud.aux = aux
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)
    th = gd_th.ThermodynamicalQuantities(ud)
    sol = fields.CellSolField(elem.sc)
    npf = fields.NodePressureField(elem, node, ud)
    sol = mod.sol_init(sol, npf, elem, node, th, ud)
    out = {}
    for attr in ("rho", "rhou", "rhov", "rhoY"):
        out[attr] = (
            np.asarray(getattr(sol, attr)[elem.i2]).reshape(elem.iicx, elem.iicy).T
        )
    out["p2_nodes"] = np.asarray(npf.p2_nodes[node.i2]).reshape(node.iicx, node.iicy).T
    return out


def t0_point(work, truth):
    """(rmse, spread) per field from the t=0 prior snapshots vs the truth IC."""
    res = {}
    for field, tru in truth.items():
        mems = sorted(
            glob.glob(f"{WT}/{work}/memory/model/pybella/*_0000/prior_mem*/{field}.npy")
        )
        if not mems:
            continue
        ens = np.array([np.load(m) for m in mems])
        rmse = float(np.sqrt(((ens.mean(0) - tru) ** 2).mean()))
        res[field] = (rmse, float(ens.std(0, ddof=1).mean()))
    return res


def load(path, state="after_full_step"):
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["state"] == state:
                d.setdefault(r["field"], []).append(
                    (float(r["t"]), float(r["rmse"]), float(r["spread"]))
                )
    return {k: sorted(v) for k, v in d.items()}


CASES = {
    "bubble": {
        "t0": ("rb", "CFLfixed_truth", "nedas_bubble_enda"),
        "spinup_csv": f"{WT}/nedas_diag_bubble_enda/rmse_spread.csv",  # pre-DA priors
        "curves": {
            "native EnDA": f"{NAT}/osse_final_bubble_enda/rmse_spread.csv",
            "native EnDAB": f"{NAT}/osse_final_bubble_endab/rmse_spread.csv",
            "NEDAS EnDA": f"{WT}/nedas_diag_bubble_enda/rmse_spread.csv",
            "NEDAS EnDAB": f"{WT}/nedas_diag_bubble_endab/rmse_spread.csv",
        },
    },
    "tv_momentum": {
        "t0": ("test_travelling_vortex", "truth", "nedas_tv_enda"),
        "spinup_csv": None,
        "curves": {
            "native EnDA": f"{NAT}/osse_final_tv_enda/rmse_spread.csv",
            "native EnDAB": f"{NAT}/osse_final_tv_endab/rmse_spread.csv",
            "NEDAS EnDA": f"{WT}/nedas_diag_tv_enda/rmse_spread.csv",
            "NEDAS EnDAB": f"{WT}/nedas_diag_tv_endab/rmse_spread.csv",
        },
    },
    "tv_all": {
        "t0": ("test_travelling_vortex", "truth", "nedas_tv_enda_all"),
        "spinup_csv": None,
        "curves": {
            "native EnDA": f"{NAT}/osse_final_tv_enda_all/rmse_spread.csv",
            "native EnDAB": f"{NAT}/osse_final_tv_endab_all/rmse_spread.csv",
            "NEDAS EnDA": f"{WT}/nedas_diag_tv_enda_all/rmse_spread.csv",
            "NEDAS EnDAB": f"{WT}/nedas_diag_tv_endab_all/rmse_spread.csv",
        },
    },
}
COLORS = {
    "native EnDA": "tab:blue",
    "native EnDAB": "tab:green",
    "NEDAS EnDA": "tab:orange",
    "NEDAS EnDAB": "tab:red",
}
outdir = f"{WT}/nedas_rmse_overlays"
os.makedirs(outdir, exist_ok=True)

for case, cfg in CASES.items():
    ic_key, aux, t0_work = cfg["t0"]
    t0 = t0_point(t0_work, truth_ic(ic_key, aux))
    data = {lbl: load(p) for lbl, p in cfg["curves"].items() if os.path.exists(p)}
    spinup = load(cfg["spinup_csv"], state="before_da") if cfg["spinup_csv"] else {}
    fields = sorted({f for d in data.values() for f in d})
    fig, axes = plt.subplots(1, len(fields), figsize=(4.2 * len(fields), 3.6))
    axes = np.atleast_1d(axes)
    for ax, field in zip(axes, fields):
        for label, d in data.items():
            if field not in d:
                continue
            rows = d[field]
            if field in t0 and t0[field][0] > 1e-12:
                rows = [(0.0, *t0[field])] + rows
            ts = [r[0] for r in rows]
            ax.plot(
                ts,
                [r[1] for r in rows],
                "-o",
                ms=3,
                color=COLORS[label],
                label=f"{label} RMSE",
            )
            ax.plot(
                ts,
                [r[2] for r in rows],
                "--",
                lw=1,
                color=COLORS[label],
                label=f"{label} spread",
            )
        if field in spinup:
            pre = [(r) for r in spinup[field] if r[0] < 0.5]
            if field in t0 and t0[field][0] > 1e-12:
                pre = [(0.0, *t0[field])] + pre
            ax.plot(
                [r[0] for r in pre],
                [r[1] for r in pre],
                ":",
                lw=1.4,
                color="gray",
                label="spinup fcst (NEDAS)",
            )
        if field in t0 and t0[field][0] > 1e-12:
            ax.plot(
                [0.0], [t0[field][0]], "k*", ms=9, zorder=5, label="initial error (t=0)"
            )
        ax.set_yscale("log")
        ax.set_xlabel("t")
        ax.set_title(field)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("analysis RMSE (solid) / spread (dashed)")
    axes[-1].legend(fontsize=6, ncol=1, loc="best")
    fig.suptitle(
        f"{case}: native LETKF vs NEDAS ETKF (analysis state, incl. t=0 initial error)",
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = f"{outdir}/rmse_overlay_{case}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out, "| t0:", {k: f"{v[0]:.2e}" for k, v in t0.items()})
