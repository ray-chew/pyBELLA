"""Time series + reflection diagnostics for study.py runs.

    python analyze.py OUT_PNG_PREFIX label1=run1.npz label2=run2.npz ...

1. ``<prefix>_series.png``: rel-L2 (w, u') vs the strict oracle in the
   windows surf 0-1 km / 1-4.5 km / full 0-12 km, and the drag ratio,
   against t U / a.
2. ``<prefix>_reflection.png`` + stdout: per-wavenumber downward/upward
   amplitude ratio |B/A| from a least-squares fit of the density-scaled
   eta-level spectrum, w_hat(k, z) sqrt(rho/rho0) = A e^{i m z} + B e^{-i m z}
   (m = sqrt(l^2 - k^2)), over 0-12 km, for every saved snapshot — the
   direct measure of top/sponge reflection. The oracle's own |B/A| is
   printed as the floor (non-Boussinesq structure is not exactly e^{imz}).
"""

import ast
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pybella.tests import linear_mountain_oracle as lo  # noqa: E402

U, N = 10.0, 0.01


def series(d, key):
    keys = list(d["series_keys"])
    return d["series"][:, keys.index(key)]


def reflection(w, zeta, x, a_hill, zmax=12000.0):
    """|B/A| per k (k < l) on eta rows (flat away from the hill)."""
    bg = lo.Background(N_const=N)
    rows = (zeta > 300.0) & (zeta <= zmax)
    zz = zeta[rows]
    sc = np.sqrt(bg.rho(zz) / bg.rho(0.0))
    nx = x.size
    k = 2 * np.pi * np.fft.fftfreq(nx, d=x[1] - x[0])
    W = np.fft.fft(w[:, rows] * sc[None, :], axis=0)
    l = N / U
    out = []
    for i in range(1, nx // 2):
        if k[i] >= 0.9 * l:
            break
        m = np.sqrt(l**2 - k[i] ** 2)
        M = np.stack([np.exp(1j * m * zz), np.exp(-1j * m * zz)], axis=1)
        coef, *_ = np.linalg.lstsq(M, W[i], rcond=None)
        out.append((k[i], abs(coef[1]) / max(abs(coef[0]), 1e-30), abs(coef[0])))
    return np.array(out)


def main():
    prefix = sys.argv[1]
    runs = [a.split("=", 1) for a in sys.argv[2:]]
    data = {lab: np.load(fn) for lab, fn in runs}

    fig, axs = plt.subplots(2, 4, figsize=(18, 7), sharex=True)
    for lab, d in data.items():
        args = ast.literal_eval(d["args"].item())
        a_hill = 1000.0 if args["case"] == "narrow" else 10000.0
        t = series(d, "t") * U / a_hill
        for j, win in enumerate(
            ("surf_0-1km", "old_1-4.5km", "mid_4.5-12km", "full_0-12km")
        ):
            axs[0, j].plot(t, series(d, f"{win}/w"), label=lab)
            axs[1, j].plot(t, series(d, f"{win}/u"), label=lab)
            axs[0, j].set_title(f"w rel-L2  {win}")
            axs[1, j].set_title(f"u' rel-L2  {win}")
    for ax in axs.ravel():
        ax.set_ylim(0, 1.2)
        ax.grid(alpha=0.3)
    for ax in axs[1]:
        ax.set_xlabel("t U / a")
    axs[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "_series.png", dpi=110)
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    for lab, d in data.items():
        args = ast.literal_eval(d["args"].item())
        a_hill = 1000.0 if args["case"] == "narrow" else 10000.0
        t = series(d, "t") * U / a_hill
        ax.plot(t, series(d, "drag_ratio_1-9km"), label=lab)
    ax.axhline(1.0, color="k", lw=0.5)
    ax.set_xlabel("t U / a")
    ax.set_ylabel("flux(1-9 km) / D_oracle")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "_drag.png", dpi=110)
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    for lab, d in data.items():
        args = ast.literal_eval(d["args"].item())
        a_hill = 1000.0 if args["case"] == "narrow" else 10000.0
        x, zeta = d["x"], d["zeta"]
        r0 = reflection(d["ref_w"], zeta, x, a_hill)
        print(
            f"{lab}: oracle floor |B/A| (k-weighted) "
            f"{np.sum(r0[:, 1] * r0[:, 2]) / np.sum(r0[:, 2]):.3f}"
        )
        for i, ts in enumerate(d["snap_times"]):
            r = reflection(d[f"sim{i}_w"], zeta, x, a_hill)
            wavg = np.sum(r[:, 1] * r[:, 2]) / np.sum(r[:, 2])
            print(
                f"   t={ts:7.0f}  amplitude-weighted |B/A| = {wavg:.3f}   "
                f"per k(<=6): " + " ".join(f"{v:.2f}" for v in r[:6, 1])
            )
            if i == len(d["snap_times"]) - 1:
                ax.semilogx(r[:, 0], r[:, 1], "o-", label=f"{lab} t={ts:.0f}")
        ax.semilogx(r0[:, 0], r0[:, 1], "k--", lw=0.8)
    ax.set_xlabel("k [1/m]")
    ax.set_ylabel("|B / A| (down / up)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "_reflection.png", dpi=110)
    plt.close(fig)
    print("wrote", prefix + "_{series,drag,reflection}.png")


if __name__ == "__main__":
    main()
