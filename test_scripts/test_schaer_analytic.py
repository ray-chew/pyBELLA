"""Schär ridge: linear-oracle gates + the Gal-Chen-vs-SLEVE discriminator.

Runs the ``test_schaer_ridge`` configuration in-process on a reduced grid
(+-32 km, 128x48, dx = 500 m) to quasi-steady state (480 x 12.5 s = 6000 s,
t U / a = 12) twice — once with the case's SLEVE transform, once forced to
Gal-Chen — and asserts:

1. SLEVE vs the linear FFT oracle (``tests/schaer_linear_analytic.py``):
   field, drag and flux-constancy gates. Calibration: w 0.312, u' 0.307,
   drag_ratio 1.104, flux_constancy 0.043. (At 800 steps the
   periodic-domain mean-flow deceleration drags drag_ratio to ~1.2 — the
   metrics are quoted at the 480-step quasi-steady window on purpose.)
2. The discriminator: the small-scale spectral fraction of w at 4-9 km —
   where the lambda = 4 km response is evanescent-dead in the true
   solution — must collapse under SLEVE relative to Gal-Chen
   (calibration: E_ss 0.004 vs 0.070, a factor ~17) and w must agree with
   linear theory at least as well. Relative assertions: calibration drift
   cannot silently invert the conclusion.

The oracle itself is self-tested here against Smith's closed-form Agnesi
drag and the evanescent-ridge zero-drag property (no simulation needed).
"""

import matplotlib

matplotlib.use("Agg")

import os

import matplotlib.pyplot as plt
import numpy as np

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import terrain
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import agnesi_smith_analytic as smith
from pybella.tests import schaer_linear_analytic as lin
from pybella.tests import test_schaer_ridge as case
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState


class _StubWriter:
    def write(self, *args, **kwargs):
        pass

    def populate(self, *args, **kwargs):
        pass

    def write_all(self, *args, **kwargs):
        pass


def run_to_steady_state(transform="sleve", steps=480):
    """Reduced oracle config: +-32 km, 128x48 (dx 500 m), sponge above
    ~11.4 km, t U / a = 12. ~10 s per run."""
    ud = user_data.UserDataInit(**vars(case.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.xmin, ud.xmax = -32000.0 / ud.h_ref, 32000.0 / ud.h_ref
    ud.inx = 128 + 1
    ud.iny = 48 + 1
    ud.inbcy = 20
    ud.stepmax = steps
    ud.tout = [1e6]
    ud.diag = False
    if transform == "galchen":
        ud.vertical_transform = terrain.GalChenTransform()
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = case.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    mem = time_update.do(mem, ud, tout=ud.tout[0], debug_writer=_StubWriter())
    return mem, ud


# --- oracle self-tests (no simulation) ---------------------------------------


def test_oracle_reproduces_smith_drag():
    U, N, h0, a, rho0 = 10.0, 0.01, 100.0, 10000.0, 1.2
    L, nx = 400e3, 2048
    x = (np.arange(nx) + 0.5) * (L / nx) - L / 2
    h = h0 * a**2 / (x**2 + a**2)
    D = lin.analytic_drag(x, h, U, N, rho0)
    D_smith = 0.25 * np.pi * rho0 * N * U * h0**2
    assert abs(D / D_smith - 1.0) < 0.02


def test_oracle_evanescent_ridge_has_no_drag():
    U, N = 10.0, 0.01
    L, nx = 400e3, 2048
    x = (np.arange(nx) + 0.5) * (L / nx) - L / 2
    h = 50.0 * np.cos(2 * np.pi * x / 4000.0) * np.exp(-((x / 20000.0) ** 2))
    D = lin.analytic_drag(x, h, U, N, 1.2)
    D_scale = 0.25 * np.pi * 1.2 * N * U * 100.0**2
    assert abs(D) < 1e-6 * D_scale


def test_oracle_fields_match_smith_for_agnesi():
    """Analytic-vs-analytic: the FFT solution against Smith's closed form
    (hydrostatic approximation) for the Agnesi profile — agreement to the
    size of the nonhydrostatic correction at N a / U = 10."""
    U, N, h0, a = 10.0, 0.01, 100.0, 10000.0
    L, nx = 400e3, 1024
    x = (np.arange(nx) + 0.5) * (L / nx) - L / 2
    h = h0 * a**2 / (x**2 + a**2)
    zlev = np.linspace(500.0, 9000.0, 30)
    z = np.broadcast_to(zlev, (nx, zlev.size)).copy()
    w_f, _ = lin.linear_fields(x, z, h, U, N)
    _, w_s, _ = smith.smith_fields(x, zlev, {"U": U, "N": N, "h0": h0, "a": a})
    assert np.linalg.norm(w_f - w_s) / np.linalg.norm(w_s) < 0.2


# --- simulation gates + discriminator ----------------------------------------


def _plot_comparison(metrics_by_tr, fields_by_tr, ud):
    out_dir = os.path.join("outputs", "schaer_discriminator")
    os.makedirs(out_dir, exist_ok=True)
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.2), sharey=True)
    (x, z, w_ref), labels = fields_by_tr["ref"], ["sleve", "galchen"]
    for ax, key in zip(axs[:2], labels):
        xx, zz, w = fields_by_tr[key]
        lim = np.abs(w_ref).max()
        ax.pcolormesh(
            xx / 1e3, zz / 1e3, w, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto"
        )
        ax.set_title(f"w [{key}]  E_ss={metrics_by_tr[key]['E_ss']:.3f}")
        ax.set_xlabel("x [km]")
    lim = np.abs(w_ref).max()
    axs[2].pcolormesh(
        x / 1e3, z / 1e3, w_ref, cmap="RdBu_r", vmin=-lim, vmax=lim, shading="auto"
    )
    axs[2].set_title("w [linear oracle]")
    axs[2].set_xlabel("x [km]")
    axs[0].set_ylabel("z [km]")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "schaer_w_discriminator.png"), dpi=130)
    plt.close(fig)


def test_schaer_discriminator():
    metrics = {}
    plot_fields = {}
    for tr in ("sleve", "galchen"):
        mem, ud = run_to_steady_state(tr)
        m, _ = lin.compare(mem, ud, z_lo_SI=1000.0, z_hi_SI=6000.0)
        metrics[tr] = m
        x, z, _, w_sim, _ = smith.sim_perturbations_SI(mem, ud)
        xx = np.broadcast_to(x.reshape(-1, 1), z.shape)
        plot_fields[tr] = (xx, z, w_sim)

    h_SI = ud.orography(x / ud.h_ref, 0.0) * ud.h_ref
    w_ref, _ = lin.linear_fields(x, z, h_SI, ud.U0, ud.NN)
    plot_fields["ref"] = (xx, z, w_ref)
    _plot_comparison(metrics, plot_fields, ud)

    sl, gc = metrics["sleve"], metrics["galchen"]

    # 1. SLEVE vs linear oracle (calibrated 0.312 / 0.307 / 1.104 / 0.043;
    #    Nh0/U = 0.25 is weakly nonlinear — gates ~40% above calibration)
    assert sl["w"] <= 0.45, f"w rel-L2 vs oracle: {sl['w']:.3f}"
    assert sl["u"] <= 0.45, f"u' rel-L2 vs oracle: {sl['u']:.3f}"
    assert 0.90 <= sl["drag_ratio"] <= 1.25, f"drag ratio: {sl['drag_ratio']:.3f}"
    assert sl["flux_constancy"] <= 0.10, f"flux constancy: {sl['flux_constancy']:.3f}"

    # 2. the discriminator (relative — calibration drift cannot invert it):
    #    small-scale w aloft collapses under SLEVE (calibrated 0.004 vs
    #    0.070), and the wave field agrees with linear theory at least as
    #    well as under Gal-Chen
    assert sl["E_ss"] <= 0.5 * gc["E_ss"], (
        f"SLEVE small-scale fraction {sl['E_ss']:.4f} not below half of "
        f"Gal-Chen's {gc['E_ss']:.4f}"
    )
    assert sl["E_ss"] <= 0.02, f"absolute small-scale fraction: {sl['E_ss']:.4f}"
    assert (
        sl["w"] <= gc["w"] + 0.02
    ), f"SLEVE w ({sl['w']:.3f}) worse than Gal-Chen ({gc['w']:.3f})"
