"""Reproduce the Chew-Benacchio-Klein 2022 (MWR) OSSE experiments.

Runs, per case, the five-run recipe of the paper-era queue driver
(archive/localdab:RKLM_Python/queue_run.py), extended with the ETPF variants:

1. obs    N=1  observation-generating run (seeded truth IC)
2. truth  N=1  truth run for the RMSE reference (initial blending on)
3. noda   N=K  ensemble forecast, no assimilation
4. enda   N=K  LETKF (rloc), momentum observations, no blending
5. endab  N=K  LETKF (rloc) + one blended step after each assimilation
6. etpf   N=K  ETPF, no blending
7. etpfb  N=K  ETPF + blending

Observation files are regenerated from seeds, never committed: the obs run's
output HDF5 *is* the observation file (labels ``<attr>_ensemble_mem=0_<t>_
after_full_step`` are exactly what ``da_params.init.load_obs`` expects).

Usage (from the repository root):

    python run_scripts/osse_mwr2022.py tv
    python run_scripts/osse_mwr2022.py bubble --members 10 --runs obs truth enda
"""

import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from driver import run_params  # noqa: E402

CASES = {
    "tv": {
        "ic": "test_travelling_vortex",
        # thesis Table 6.2: t_first = 0.25, dt_obs = 0.25, t_end = 3.0
        "tout": [round(float(t), 3) for t in np.arange(0.25, 3.25, 0.25)],
        "ud_extra": {},
        # TV sol_init keys the (identical) truth/obs IC on 'obs'/'truth' in aux
        "obs_aux": "obs",
        "truth_aux": "truth",
    },
    "bubble": {
        # the production paper case (t_ref = 1000 s, native 160x80, seeded
        # delth machinery); tests/test_blending_warm_bubble is only a 31-step
        # blending smoke and is unstable at the paper grid/times
        "ic": "rb",
        # thesis Table 6.2: t_first = 500 s, dt_obs = 50 s, t_end = 1000 s,
        # nondimensionalised by t_ref
        "tout": [round(float(t), 3) for t in np.arange(0.5, 1.05, 0.05)],
        "ud_extra": {},
        # the bubble sol_init keys the truth IC on 'truth' in aux only
        "obs_aux": "truth_obs",
        "truth_aux": "truth",
        # rb keys behaviour on the aux substring 'CFLfixed' (2-step dt pin).
        # Do NOT include 'imbal': it triggers the initial *hydrostatic*
        # conversion (hydrostatic-blending experiments), which is unstable for
        # this nonhydrostatic case; the paper's initial blending is the
        # pseudo-incompressible one driven by initial_blending=True alone.
        "aux_base": "CFLfixed_",
        "out_glob": "./outputs/output_rising_bubble/*ensemble=1*_%s*.h5",
    },
}

# thesis Table 6.2: momentum-only observations
OBS_ATTRS = ["rhou", "rhov"]

ALL_RUNS = ["obs", "truth", "noda", "enda", "endab", "etpf", "etpfb"]


def find_obs_file(case_cfg, grid=None):
    default = "./outputs/test_*%s*/*ensemble=1*_%%s*.h5" % (
        case_cfg["ic"].replace("test_", "")
    )
    pattern = case_cfg.get("out_glob", default) % case_cfg["obs_aux"]
    if grid is not None:
        # pin the grid so e.g. a low-resolution obs file cannot shadow the
        # paper-resolution one (fn_gen writes _<Nx>_<Ny>_ into the suffix)
        pattern = pattern.replace("*ensemble=1*", "*ensemble=1_%i_%i*" % grid)
    matches = sorted(glob.glob(pattern))
    assert matches, "no observation file found for pattern %s -- run 'obs' first" % (
        pattern
    )
    return matches[-1]


def queue(ic, N, ud, dap):
    rp = run_params()
    rp.tc = ic
    rp.N = N
    rp.ud = json.dumps(ud)
    rp.dap = json.dumps(dap)
    print(">>> pybella -ic %s -N %i queue -w '%s' '%s'" % (ic, N, rp.ud, rp.dap))
    rp.queue_run()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("case", choices=list(CASES))
    ap.add_argument("--members", type=int, default=10, help="ensemble size K")
    ap.add_argument("--runs", nargs="*", default=ALL_RUNS, choices=ALL_RUNS)
    ap.add_argument(
        "--ud",
        type=json.loads,
        default={},
        help="extra ud overrides merged into every run, e.g. "
        '\'{"inx": 81, "iny": 41}\' for a cheap shakedown',
    )
    args = ap.parse_args()

    cfg = CASES[args.case]
    ic, K = cfg["ic"], args.members

    ud_common = {
        "diag": False,  # CompareSol would gate every member against the N=1 target
        "autogen_fn": True,
        "stepmax": 100000,
        "tout": cfg["tout"],
    }
    ud_common.update(cfg["ud_extra"])
    ud_common.update(args.ud)

    grid = None
    if "inx" in ud_common and "iny" in ud_common:
        grid = (ud_common["inx"] - 1, ud_common["iny"] - 1)

    aux_base = cfg.get("aux_base", "")

    def ud_for(aux, blend=False):
        ud = dict(ud_common, aux=aux_base + aux, initial_blending=True)
        if blend:
            # one blended (pseudo-incompressible) step after each assimilation
            ud["continuous_blending"] = True
        return ud

    dap_noda = {"da_times": []}

    def dap_for(da_type):
        # resolved lazily: the obs file may be produced by the 'obs' run of
        # this same invocation
        return {
            "da_times": cfg["tout"],
            "obs_attrs": OBS_ATTRS,
            "obs_path": find_obs_file(cfg, grid=grid),
            "da_type": da_type,
        }

    if "obs" in args.runs:
        # no initial blending for the observation run (archive queue_run.py)
        queue(ic, 1, dict(ud_common, aux=aux_base + cfg["obs_aux"]), dap_noda)
    if "truth" in args.runs:
        queue(ic, 1, ud_for(cfg["truth_aux"]), dap_noda)
    if "noda" in args.runs:
        queue(ic, K, ud_for("noda"), dap_noda)
    if "enda" in args.runs:
        queue(ic, K, ud_for("wda"), dap_for("rloc"))
    if "endab" in args.runs:
        queue(ic, K, ud_for("wda", blend=True), dap_for("rloc"))
    if "etpf" in args.runs:
        queue(ic, K, ud_for("wda_etpf"), dap_for("etpf"))
    if "etpfb" in args.runs:
        queue(ic, K, ud_for("wda_etpf", blend=True), dap_for("etpf"))


if __name__ == "__main__":
    main()
