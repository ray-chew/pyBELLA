"""Physics oracle: the igw_baldauf_brdar case vs the linear analytic reference.

Runs the Baldauf & Brdar internal-gravity-wave case in-process (regression
configuration: dt = 500 s, dx ~ 20 km, t_end = 15500 s, f-plane Coriolis) and
compares the final x-mean-free wave fields against the numerically-exact
linear reference (``pybella.tests.baldauf_brdar_analytic``). Unlike the
golden-master comparisons (which catch *change*), this catches *wrongness*:
a flipped Coriolis or buoyancy sign, a broken dispersion relation, or wrong
wave amplitudes blow these bounds by an order of magnitude (a pure sign
error alone gives rel L2 ~ 1.4-2.0).

Measured rel-L2 (2026-06-10, regression config, AFTER the 2D out-of-plane
Coriolis fixes in BOTH explicit_euler.do_forward_step and the implicit
correction's w-row): u 0.285, vo 0.045, w 0.556, p 0.670, rho 0.346.
Gates are >=1.4x above. The residual was shown
to decompose into (validated by refinement studies, see
dev_notes/regression_harness.md):

- temporal:  O((omega dt)^2) phase error — at dt = 125 s:
  u 0.133, vo 0.028, w 0.272, rho 0.221;
- spatial:   under-resolved w/rho structure at dx = 20 km — at dx ~ 10 km
  (f = 0, dt = 125 s): u 0.080, w 0.166, rho 0.197;
- p' is split between rhoY and p2 internally at these amplitudes, so p is
  only loosely gated.

History: this oracle originally measured vo pinned at 0.44 rel-L2
independent of dt with a sim/ref amplitude ratio ~0.6 — the fingerprint of
the (since fixed) 2D defect where the rhow row of the explicit forward step
sat behind ``if ndim == 3`` and the out-of-plane momentum received only the
implicit half of the Coriolis rotation.

Comparator internals validated separately: energy drift ~1e-13 (neutral
discretisation, exact-in-time eigenpropagation), z-refinement converged
(refine 2 vs 4 identical to 3 digits), t = 0 round-trip ~3e-4.

Also writes ref/sim/diff snapshot PNGs for u' and w' into the case's output
directory so the wave physics can be verified by eye alongside the other
regression images.
"""

import numpy as np

from pybella.tests import baldauf_brdar_analytic as bb

GATES = {
    "u": 0.40,
    "vo": 0.10,
    "w": 0.80,
    "p": 1.00,  # loose: pressure decomposition is scheme-internal
    "rho": 0.50,
}


def _demean(q):
    return q - q.mean(axis=0, keepdims=True)


def test_igw_matches_linear_analytic():
    ud, ic, end, zc, t_end = bb.run_sim()

    par = bb.IGWParams(ud)
    L = (ud.xmax - ud.xmin) * ud.h_ref
    ref, diag = bb.evolve_linear(ic, L, zc, t_end, par)

    # comparator-internal exactness
    assert diag["energy_drift"] < 1e-9

    failures = []
    for key, gate in GATES.items():
        s, r = _demean(end[key]), ref[key]
        rel = np.linalg.norm(s - r) / np.linalg.norm(r)
        if rel > gate:
            failures.append(f"{key}: rel L2 {rel:.3f} > gate {gate}")

    # amplitude sanity on the primary wave field
    amp_ratio = np.abs(_demean(end["u"])).max() / np.abs(ref["u"]).max()
    if not (0.7 < amp_ratio < 1.3):
        failures.append(f"u amplitude ratio sim/ref {amp_ratio:.2f} outside (0.7, 1.3)")

    _write_pngs(end, ref, zc, L)

    assert not failures, "; ".join(failures)


def _write_pngs(end, ref, zc, L):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = "./outputs/test_igw_baldauf_brdar/"
    nx = end["u"].shape[0]
    x = (np.arange(nx) + 0.5) * (L / nx) / 1000.0 - L / 2000.0  # km
    zk = zc / 1000.0

    for key, label in (("u", "u' [m/s]"), ("w", "w' [m/s]")):
        s, r = _demean(end[key]), ref[key]
        fig, axs = plt.subplots(1, 3, figsize=(14, 3), sharey=True)
        for ax, (arr, title) in zip(
            axs, ((r, "linear analytic"), (s, "pyBELLA"), (s - r, "diff"))
        ):
            pc = ax.pcolormesh(x, zk, arr.T, shading="auto", cmap="RdBu_r")
            ax.set_title(f"{label} — {title}")
            ax.set_xlabel("x [km]")
            fig.colorbar(pc, ax=ax)
        axs[0].set_ylabel("z [km]")
        fig.tight_layout()
        fig.savefig(outdir + f"analytic_{key}.png", dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    test_igw_matches_linear_analytic()
