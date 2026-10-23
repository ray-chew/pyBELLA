"""Mountain-wave model-vs-strict-linear-oracle study driver.

Runs a terrain case in-process (native 2D), and at regular snapshots
compares the model against ``tests/linear_mountain_oracle.py`` (the oracle
is time-independent, so it is evaluated once). Writes an ``.npz`` with the
metric time series, field snapshots and the oracle fields.

    python run_scripts/strict_linear_oracle/study.py --case agnesi \
        --nx 96 --nz 48 --dt 50 --tend 40000 --out outputs/strict_oracle/a.npz

Cases: ``agnesi`` (test_agnesi_hydrostatic: a = 10 km, N a/U = 10),
``narrow`` (same, a = 1 km, N a/U = 1, +-40 km, 20 km top, sponge > 10 km),
``schaer`` (test_schaer_ridge, SLEVE). ``--sponge-rate-fixed`` (default on)
rescales the per-step Rayleigh coefficients so the damping *rate* is that
of dt = 50 s, independent of dt.
"""

import argparse
import os
import time

import numpy as np

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import terrain
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import agnesi_smith_analytic as smith
from pybella.tests import linear_mountain_oracle as lo
from pybella.tests import test_agnesi_hydrostatic as agnesi_case
from pybella.tests import test_schaer_ridge as schaer_case
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState


class _StubWriter:
    def write(self, *a, **k):
        pass

    def populate(self, *a, **k):
        pass

    def write_all(self, *a, **k):
        pass


WINDOWS = {
    "surf_0-1km": (0.0, 1000.0),
    "old_1-4.5km": (1000.0, 4500.0),
    "mid_4.5-12km": (4500.0, 12000.0),
    "full_0-12km": (0.0, 12000.0),
}


def build_ud(args):
    case = schaer_case if args.case == "schaer" else agnesi_case
    # hill parameters must be set on the UserData INSTANCE: ud.orography is
    # its bound method and reads self.hill_height / self.hill_width
    inst = case.UserData()
    if args.case == "narrow":
        inst.hill_width = 1000.0
    if args.h0 is not None:
        inst.hill_height = args.h0
    ud = user_data.UserDataInit(**vars(inst))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    if args.case == "narrow":
        ud.xmin, ud.xmax = -40e3 / ud.h_ref, 40e3 / ud.h_ref
        ud.ymax = 20e3 / ud.h_ref
    if args.case == "schaer":
        ud.xmin, ud.xmax = -32e3 / ud.h_ref, 32e3 / ud.h_ref
    if args.ymax is not None:
        ud.ymax = args.ymax / ud.h_ref
    if args.xhalf is not None:
        ud.xmin, ud.xmax = -args.xhalf / ud.h_ref, args.xhalf / ud.h_ref
    assert np.isclose(ud.orography(0.0, 0.0) * ud.h_ref, ud.hill_height)
    if args.galchen:
        ud.vertical_transform = terrain.GalChenTransform()
    ud.inx = args.nx + 1
    ud.iny = args.nz + 1
    ud.inz = 1
    dy = (ud.ymax - ud.ymin) * ud.h_ref / args.nz
    ud.inbcy = int(round((ud.ymax * ud.h_ref - args.sponge_base) / dy))
    ud.dtfixed = ud.dtfixed0 = args.dt / ud.t_ref
    ud.stepmax = 10**9
    ud.tout = [1e9]
    ud.diag = False
    if args.hydro:
        ud.is_nonhydrostatic = 0
    return ud, case


def rescale_sponge(ud, dt, scale=1.0):
    """Per-step factor (1 + tau) -> (1 + tau)^(scale dt / 50 s): the damping
    RATE of the code's sponge at dt = 50 s, times ``scale``."""
    r = scale * dt / 50.0
    ud.tcy = (1.0 + ud.tcy) ** r - 1.0
    ud.tny = (1.0 + ud.tny) ** r - 1.0


def sin2_sponge(ud, elem, node, dt, nu_max, z_base):
    """Klemp-Lilly-type sponge: rate nu(z) = nu_max sin^2(pi/2 s),
    s = (eta - z_base) / (top - z_base), on the eta coordinate; per-step
    factor exp(-nu dt)."""
    top = ud.ymax * ud.h_ref
    for attr, coords in (("tcy", elem.y), ("tny", node.y)):
        zz = coords * ud.h_ref
        sfrac = np.clip((zz - z_base) / (top - z_base), 0.0, 1.0)
        nu = nu_max * np.sin(0.5 * np.pi * sfrac) ** 2
        setattr(ud, attr, np.exp(-nu * dt) - 1.0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--case", default="agnesi", choices=("agnesi", "narrow", "schaer"))
    p.add_argument("--nx", type=int, default=96)
    p.add_argument("--nz", type=int, default=48)
    p.add_argument("--dt", type=float, default=50.0)
    p.add_argument("--tend", type=float, default=24000.0)
    p.add_argument("--snap", type=float, default=2000.0, help="metric interval [s]")
    p.add_argument("--save-fields", type=str, default="", help="times [s], comma sep")
    p.add_argument("--h0", type=float, default=None)
    p.add_argument("--sponge-base", type=float, default=12000.0)
    p.add_argument("--hydro", action="store_true", help="is_nonhydrostatic = 0")
    p.add_argument("--galchen", action="store_true")
    p.add_argument("--no-sponge-rescale", action="store_true")
    p.add_argument("--sponge-scale", type=float, default=1.0)
    p.add_argument(
        "--sponge-numax",
        type=float,
        default=0.0,
        help="> 0: replace the code sponge by a sin^2 rate profile [1/s]",
    )
    p.add_argument("--ymax", type=float, default=None, help="domain top [m]")
    p.add_argument("--xhalf", type=float, default=None, help="domain half-width [m]")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    ud, case = build_ud(args)
    zs = args.sponge_base
    WINDOWS["mid_4.5-12km"] = (4500.0, zs)
    WINDOWS["full_0-12km"] = (0.0, zs)  # keys kept; upper edge = sponge base
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = case.sol_init(sol, npf, elem, node, th, ud)
    if args.sponge_numax > 0:
        sin2_sponge(ud, elem, node, args.dt, args.sponge_numax, args.sponge_base)
    elif not args.no_sponge_rescale:
        rescale_sponge(ud, args.dt, args.sponge_scale)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)

    sf0 = lo.sim_fields(mem, ud)
    dx_SI = elem.dx * ud.h_ref
    t0 = time.time()
    regime_main = dict(regime="compressible", hydrostatic=bool(args.hydro))
    ref = lo.oracle_for_sim(sf0, ud, **regime_main)
    alt = {
        "compNH": lo.oracle_for_sim(sf0, ud, "compressible", False),
        "compH": lo.oracle_for_sim(sf0, ud, "compressible", True),
        "bousH": lo.oracle_for_sim(sf0, ud, "boussinesq", True),
        "bousNH": lo.oracle_for_sim(sf0, ud, "boussinesq", False),
    }
    bg = lo.Background.from_ud(ud)
    h = ud.orography(sf0["x"] / ud.h_ref, 0.0) * ud.h_ref
    D_ref = lo.drag_profile(sf0["x"], [0.0], h, ud.U0, bg, **regime_main)[0]
    print(
        f"oracle built in {time.time() - t0:.1f}s; D_ref = {D_ref:.2f} N/m", flush=True
    )

    save_t = [float(s) for s in args.save_fields.split(",") if s]
    series = []
    snaps = {}
    t_next = args.snap
    tw = time.time()
    while mem.time.t * ud.t_ref < args.tend - 1e-6:
        t_target = min(t_next, args.tend)
        mem = time_update.do(
            mem, ud, tout=t_target / ud.t_ref, debug_writer=_StubWriter()
        )
        t_SI = mem.time.t * ud.t_ref
        sf = lo.sim_fields(mem, ud)
        if not np.all(np.isfinite(sf["w"])):
            print("NaN at t =", t_SI, flush=True)
            break
        row = {"t": t_SI, "step": mem.time.step}
        m = lo.compare(sf, ref, WINDOWS, dx_SI)
        for wn in WINDOWS:
            for q, v in m[wn].items():
                row[f"{wn}/{q}"] = v
        for an, af in alt.items():
            ma = lo.compare(
                sf,
                af,
                {"full": WINDOWS["full_0-12km"], "old": WINDOWS["old_1-4.5km"]},
                dx_SI,
            )
            row[f"{an}/full/w"] = ma["full"]["w"]
            row[f"{an}/old/w"] = ma["old"]["w"]
            row[f"{an}/full/u"] = ma["full"]["u"]
        zeta = sf["zeta"]
        band = (zeta >= 1000.0) & (zeta <= 9000.0)
        row["drag_ratio_1-9km"] = float(np.mean(m["flux_sim"][band]) / D_ref)
        row["flux_constancy"] = float(
            np.std(m["flux_sim"][band]) / abs(np.mean(m["flux_sim"][band]))
        )
        surf = zeta < 1000.0
        row["drag_ratio_0-1km"] = float(np.mean(m["flux_sim"][surf]) / D_ref)
        row["umean_surf"] = float(sf["u"][:, 0].mean())
        if args.case != "schaer":
            om, _ = smith.compare(mem, ud, z_lo_SI=1000.0, z_hi_SI=4500.0)
            row["oldSmith/w"], row["oldSmith/u"] = om["w"], om["u"]
            row["oldSmith/drag"] = om["drag_ratio"]
        series.append(row)
        print(
            f"t={t_SI:7.0f} tU/a={t_SI * ud.U0 / ud.hill_width:5.1f} "
            f"full w {row['full_0-12km/w']:.3f} u {row['full_0-12km/u']:.3f} | "
            f"surf w {row['surf_0-1km/w']:.3f} u {row['surf_0-1km/u']:.3f} th "
            f"{row['surf_0-1km/theta']:.3f} pi {row['surf_0-1km/pi']:.3f} | "
            f"old w {row['old_1-4.5km/w']:.3f} | drag {row['drag_ratio_1-9km']:.3f} "
            f"| wall {time.time() - tw:.0f}s",
            flush=True,
        )
        for ts in save_t:
            if abs(t_SI - ts) < 0.5 * args.snap and ts not in snaps:
                snaps[ts] = {q: sf[q].copy() for q in ("u", "w", "theta", "pi")}
        if t_next <= mem.time.t * ud.t_ref + 1e-6:
            t_next += args.snap
    # always keep the final state
    snaps[float(mem.time.t * ud.t_ref)] = {
        q: sf[q].copy() for q in ("u", "w", "theta", "pi")
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    keys = sorted(series[0].keys())
    payload = {
        "series_keys": np.array(keys),
        "series": np.array([[r.get(k, np.nan) for k in keys] for r in series]),
        "x": sf0["x"],
        "z": sf0["z"],
        "zeta": sf0["zeta"],
        "xn": sf0["xn"],
        "zn": sf0["zn"],
        "D_ref": D_ref,
        "args": str(vars(args)),
        "snap_times": np.array(sorted(snaps)),
    }
    for q in ("u", "w", "theta", "pi"):
        payload[f"ref_{q}"] = ref[q] if q != "pi" else ref["pi_nodes"]
        for an, af in alt.items():
            payload[f"{an}_{q}"] = af[q] if q != "pi" else af["pi_nodes"]
        for i, ts in enumerate(sorted(snaps)):
            payload[f"sim{i}_{q}"] = snaps[ts][q]
    np.savez_compressed(args.out, **payload)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
