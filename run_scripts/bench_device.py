"""Backend benchmark: step wall time across backends and resolutions.

Usage (CPU box or GPU node — device selection is JAX's default):

    python run_scripts/bench_device.py --backend jax-device --res 256 --steps 50
    python run_scripts/bench_device.py --backend numpy --res 256 --steps 50
    python run_scripts/bench_device.py --backend jax-device --res 128 --case igw3d

Cases: a scaled travelling vortex (2D, periodic, elliptic-dominated) and a
scaled 3D internal-wave configuration (gravity + WALL). Reports compile/
first-step time, then median and p90 wall time per step over the remaining
steps, plus steps/s. Regression grids are too small to benchmark — use
>= 256^2 / 128^3 to measure anything meaningful; on GPUs also verify with
`jax.profiler.trace` that recurring transfers are only the dt scalar
(plus per-step writer pulls if output_timesteps is on — disable it for
benchmarks; this script does).
"""

import argparse
import sys
import time

import numpy as np


def build_vortex(res):
    from pybella.utils import user_data
    from pybella.tests import test_travelling_vortex_3d_coriolis as tv3d

    d = vars(tv3d.UserData())
    d["inx"] = res + 1
    d["iny"] = res + 1
    d["inz"] = 1
    ud = user_data.UserDataInit(**d)
    ud.coriolis_strength = np.zeros(3)
    ud.gravity_strength = np.zeros(3)
    ud.nonhydrostasy = 1.0
    ud.is_compressible = 1
    ud.compressibility = 1.0
    ud.output_timesteps = False
    ud.diag = False
    ud.aux = ""
    return ud, "xy"


def build_igw3d(res):
    from pybella.utils import user_data
    from pybella.tests import test_internal_long_wave as igw

    d = vars(igw.UserData())
    d["inx"] = res + 1
    d["iny"] = max(res // 4, 8) + 1
    d["inz"] = max(res // 4, 8) + 1
    ud = user_data.UserDataInit(**d)
    ud.output_timesteps = False
    ud.diag = False
    ud.aux = ""
    return ud, None


def build_state(ud, plane):
    from pybella.flow_solver.discretisation import grid as dis_grid
    from pybella.flow_solver.physics import hydrostatics
    from pybella.flow_solver.physics import thermodynamics
    from pybella.flow_solver.utils import cache, fields
    from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
    from pybella.utils.data_structures import ModelState

    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    hydrostatics.integrated_state(npf, elem, node, th, ud)

    rng = np.random.default_rng(11)
    shape = sol.rho.shape
    sol.rho[...] = 1.0 + 0.01 * rng.standard_normal(shape)
    sol.rhou[...] = 0.01 * rng.standard_normal(shape)
    sol.rhov[...] = 0.01 * rng.standard_normal(shape)
    sol.rhow[...] = 0.0
    sol.rhoY[...] = 1.0
    sol.rhoX[...] = 0.0

    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    if not hasattr(ud, "nonhydrostasy"):
        ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    if not hasattr(ud, "compressibility"):
        ud.compressibility = float(ud.is_compressible)
    return mem


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--backend", default="jax-device", choices=["numpy", "jax", "jax-device"]
    )
    ap.add_argument("--case", default="vortex", choices=["vortex", "igw3d"])
    ap.add_argument("--res", type=int, default=256)
    ap.add_argument("--steps", type=int, default=50)
    args = ap.parse_args()

    if args.backend != "numpy":
        import jax  # noqa: F401  (fail fast if missing)

        from pybella.backends import jax_ops  # noqa: F401  (x64)

    builder = build_vortex if args.case == "vortex" else build_igw3d
    ud, plane = builder(args.res)
    ud.backend = args.backend
    ud.stepmax = args.steps
    ud.dtfixed0 = ud.dtfixed = min(float(ud.dtfixed), 1e-3)
    mem = build_state(ud, plane)

    from pybella.flow_solver.discretisation import time_update

    class NullDebug:
        def write(self, *a, **k):
            pass

        def populate(self, *a, **k):
            pass

    # warm-up / compile: 2 steps (both Strang parities jit separately on
    # jax-device, so a 1-step warm-up would leak a compile into the timing)
    ud.stepmax = 2
    t0 = time.perf_counter()
    time_update.do(mem, ud, 1e9, None, None, NullDebug())
    compile_s = time.perf_counter() - t0

    # timed window: all steps in ONE do() call — the production pattern.
    # run_window pushes the state to device on entry and pulls it back on
    # exit, so per-call timing would charge a full host round-trip to every
    # step; window timing keeps the state device-resident throughout.
    ud.stepmax = mem.time.step + args.steps
    t0 = time.perf_counter()
    time_update.do(mem, ud, 1e9, None, None, NullDebug())
    elapsed = time.perf_counter() - t0

    per = elapsed / args.steps
    print(
        f"backend={args.backend} case={args.case} res={args.res} "
        f"steps={args.steps}\n"
        f"  warm-up (2 steps incl. compile): {compile_s:.3f} s\n"
        f"  per step: mean {per * 1e3:.2f} ms over a {args.steps}-step window, "
        f"{1.0 / per:.2f} steps/s"
    )


if __name__ == "__main__":
    sys.exit(main())
