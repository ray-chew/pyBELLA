"""Vertical profiles of (sim - oracle) at fixed x for several runs.

    python plot_profiles.py out.png q x_km snap_index zmax_km lab=run.npz ...

Distinguishes grid-scale (thickness ~ dz, numerical) from resolved
(resolution-independent, physical) near-surface structure.
"""

import ast
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def demean(q):
    return q - q.mean(axis=0, keepdims=True)


def main():
    out, q, xkm, isnap, zmax = sys.argv[1:6]
    xs = [float(v) for v in xkm.split(",")]
    isnap, zmax = int(isnap), float(zmax)
    runs = [a.split("=", 1) for a in sys.argv[6:]]
    fig, axs = plt.subplots(1, len(xs), figsize=(4.5 * len(xs), 5), squeeze=False)
    for lab, fn in runs:
        isn = isnap
        if "@" in fn:  # run.npz@snap_index overrides the default snapshot
            fn, isn = fn.split("@")
            isn = int(isn)
        d = np.load(fn)
        x, z = d["x"], d["z"]
        sim, ref = d[f"sim{isn}_{q}"], d[f"ref_{q}"]
        if q == "pi":  # nodes; drop the duplicated periodic column
            x, z, sim = d["xn"][:-1], d["zn"][:-1], sim[:-1]
        h0 = ast.literal_eval(d["args"].item())["h0"] or 100.0
        # normalised to h0 = 100 m (linear response scales with h0)
        err = (demean(sim) - demean(ref)) * 100.0 / h0
        t = d["snap_times"][isn]
        for ax, xv in zip(axs[0], xs):
            i = np.argmin(np.abs(x - xv * 1e3))
            ax.plot(err[i], z[i] / 1e3, ".-", label=f"{lab} t={t:.0f}", ms=3)
            ax.set_ylim(0, zmax)
            sel = z[i] <= zmax * 1e3
            lim = max(
                ax.get_xlim()[1] if ax.lines[:-1] else 0.0,
                1.1 * np.abs(err[i][sel]).max(),
            )
            ax.set_xlim(-lim, lim)
            ax.set_title(f"{q} (sim - oracle) x 100m/h0, x = {xv:g} km")
            ax.set_xlabel(q)
    axs[0, 0].set_ylabel("z [km]")
    axs[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
