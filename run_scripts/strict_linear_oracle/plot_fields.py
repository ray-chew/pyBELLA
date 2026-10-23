"""Field maps (sim, strict oracle, sim - oracle) from a study.py .npz.

python plot_fields.py run.npz out.png [q=w] [xlim_km=100] [zlim_km=24]
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def demean(q):
    return q - q.mean(axis=0, keepdims=True)


def main():
    fn, out = sys.argv[1], sys.argv[2]
    q = sys.argv[3] if len(sys.argv) > 3 else "w"
    xlim = float(sys.argv[4]) if len(sys.argv) > 4 else 100.0
    zlim = float(sys.argv[5]) if len(sys.argv) > 5 else 24.0
    ref_key = sys.argv[6] if len(sys.argv) > 6 else "ref"
    d = np.load(fn)
    if q == "pi":
        x, z = d["xn"][:-1], d["zn"][:-1]
    else:
        x, z = d["x"], d["z"]
    xx = np.broadcast_to(x.reshape(-1, 1), z.shape) / 1e3
    zz = z / 1e3
    ref = d[f"{ref_key}_{q}"]
    times = d["snap_times"]
    fig, axs = plt.subplots(
        len(times), 3, figsize=(15, 3.2 * len(times)), squeeze=False
    )
    inview = (np.abs(x).reshape(-1, 1) < xlim * 1e3) & (z < min(zlim, 12.0) * 1e3)
    lim = np.abs(demean(ref)[inview]).max()
    for i, t in enumerate(times):
        sim = d[f"sim{i}_{q}"]
        if q == "pi":
            sim = sim[:-1]
        for j, (fld, ttl) in enumerate(
            (
                (demean(sim), "sim"),
                (demean(ref), ref_key),
                (demean(sim) - demean(ref), "sim - oracle"),
            )
        ):
            ax = axs[i, j]
            m = ax.pcolormesh(
                xx, zz, fld, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto"
            )
            ax.set_xlim(-xlim, xlim)
            ax.set_ylim(0, zlim)
            ax.axhline(12.0, color="k", lw=0.5, ls="--")
            ax.set_title(f"{q} {ttl}  t={t:.0f}s")
        fig.colorbar(m, ax=axs[i, :].tolist(), shrink=0.8)
    for ax in axs[-1]:
        ax.set_xlabel("x [km]")
    for ax in axs[:, 0]:
        ax.set_ylabel("z [km]")
    fig.suptitle(d["args"].item()[:160], fontsize=7)
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
