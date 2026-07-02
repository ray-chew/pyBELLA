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
        "ic": "test_blending_warm_bubble",
        # thesis Table 6.2: t_first = 500 s, dt_obs = 50 s, t_end = 1000 s
        "tout": [round(float(t), 3) for t in np.arange(500.0, 1050.0, 50.0)],
        # the paper bubble is 160x80 (the case default is 64x48)
        "ud_extra": {"inx": 161, "iny": 81},
        # the bubble sol_init keys the truth IC on 'truth' in aux only
        "obs_aux": "truth_obs",
        "truth_aux": "truth",
    },
}

# thesis Table 6.2: momentum-only observations
OBS_ATTRS = ["rhou", "rhov"]

ALL_RUNS = ["obs", "truth", "noda", "enda", "endab", "etpf", "etpfb"]


def find_obs_file(case_cfg):
    pattern = "./outputs/test_*%s*/*ensemble=1*_%s*.h5" % (
        case_cfg["ic"].replace("test_", ""),
        case_cfg["obs_aux"],
    )
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

    def ud_for(aux, blend=False):
        ud = dict(ud_common, aux=aux, initial_blending=True)
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
            "obs_path": find_obs_file(cfg),
            "da_type": da_type,
        }

    if "obs" in args.runs:
        # no initial blending for the observation run (archive queue_run.py)
        queue(ic, 1, dict(ud_common, aux=cfg["obs_aux"]), dap_noda)
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
