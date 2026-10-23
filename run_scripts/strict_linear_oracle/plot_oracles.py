"""Oracle-vs-oracle figure: what the old Smith bridge misses.

For the wide (a = 10 km) and narrow (a = 1 km) Agnesi hills on the model's
periodic domains: strict compressible-NH w, and the differences
(compressible-H, Boussinesq-H x sqrt(rho0/rho) (periodic Smith + bridge),
infinite-domain Smith x sqrt(rho0/rho) (the CI comparator's reference)).

    python plot_oracles.py out.png
"""

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pybella.tests import agnesi_smith_analytic as smith  # noqa: E402
from pybella.tests import linear_mountain_oracle as lo  # noqa: E402

U, N, H0 = 10.0, 0.01, 100.0


def main():
    out = sys.argv[1]
    bg = lo.Background(N_const=N)
    cases = (
        ("wide a=10 km", 10e3, 200e3, 384, 12e3),
        ("narrow a=1 km", 1e3, 80e3, 640, 10e3),
    )
    fig, axs = plt.subplots(2, 4, figsize=(19, 7.5))
    for r, (lab, a, L, nx, ztop) in enumerate(cases):
        x = (np.arange(nx) + 0.5) * L / nx - L / 2
        h = H0 * a**2 / (x**2 + a**2)
        zl = np.linspace(0, ztop, 121)
        z = np.broadcast_to(zl, (nx, zl.size))
        cNH = lo.linear_fields(x, z, h, U, bg, "compressible", False)["w"]
        cH = lo.linear_fields(x, z, h, U, bg, "compressible", True)["w"]
        bH = lo.linear_fields(x, z, h, U, bg, "boussinesq", True)["w"]
        sc = np.sqrt(bg.rho(0.0) / bg.rho(z))
        _, wS, _ = smith.smith_fields(x, zl, {"U": U, "N": N, "h0": H0, "a": a})
        panels = (
            (cNH, "strict: compressible NH"),
            (cH - cNH, "compressible H - NH"),
            (bH * sc - cNH, "Boussinesq-H*sqrt(rho0/rho) - strict"),
            (lo.demean(wS * sc) - lo.demean(cNH), "Smith(inf)*sqrt - strict (CI ref)"),
        )
        lim = np.abs(cNH).max()
        xw = 30.0 if a > 5e3 else 15.0
        for c, (fld, ttl) in enumerate(panels):
            ax = axs[r, c]
            m = ax.pcolormesh(
                x / 1e3,
                zl / 1e3,
                fld.T,
                cmap="RdBu_r",
                vmin=-lim,
                vmax=lim,
                shading="auto",
            )
            ax.set_xlim(-xw, 2 * xw)
            ax.set_title(
                (
                    f"{lab}: {ttl}\nrel-L2 {np.linalg.norm(fld) / np.linalg.norm(cNH):.3f}"
                    if c
                    else f"{lab}: {ttl}"
                ),
                fontsize=9,
            )
            ax.set_xlabel("x [km]")
        axs[r, 0].set_ylabel("z [km]")
        fig.colorbar(m, ax=axs[r, :].tolist(), shrink=0.8, label="w [m/s]")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
