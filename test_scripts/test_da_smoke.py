"""Seeded ensemble-DA smoke test: 2 assimilation cycles on the travelling
vortex at N=4, for both LETKF (rloc) and ETPF.

Statistical assertions, not bitwise (the dask-chunked rloc analysis makes
bitwise fragile): at the final assimilation time, the analysis ensemble mean
must beat the forecast ensemble mean against the truth on the observed
momentum fields, and the analysis spread must be below the forecast spread.

The truth/observation run is regenerated here from its fixed seed (2233 via
``'obs' in ud.aux``); nothing is committed.
"""

import glob
import json
import subprocess

import h5py
import numpy as np
import pytest

MEMBERS = 4
TOUT = [0.25, 0.5]
OBS_ATTRS = ["rhou", "rhov"]
INNER = (slice(2, -2), slice(2, -2))  # igx = igy = 2


def queue_run(N, ud, dap):
    return subprocess.run(
        ["pybella", "-ic", "test_travelling_vortex", "-N", str(N), "queue", "-w"]
        + [json.dumps(ud), json.dumps(dap)],
        capture_output=True,
        text=True,
    )


def ud_payload(aux):
    return {
        "diag": False,
        "autogen_fn": True,
        "stepmax": 10000,
        "tout": TOUT,
        "aux": aux,
        "initial_blending": True,
    }


def read_field(h5, field, member, t, suffix):
    label = "%s_ensemble_mem=%i_%.3f_%s" % (field, member, t, suffix)
    return np.squeeze(h5[field][label][:])[INNER]


def rmse_and_spread(ens_h5, truth, field, t, suffix):
    members = np.array(
        [read_field(ens_h5, field, n, t, suffix) for n in range(MEMBERS)]
    )
    rmse = np.sqrt(((members.mean(axis=0) - truth) ** 2).mean())
    spread = members.std(axis=0, ddof=1).mean()
    return rmse, spread


@pytest.fixture(scope="module")
def obs_file():
    result = queue_run(
        1, dict(ud_payload("ci_obs"), initial_blending=False), {"da_times": []}
    )
    assert result.returncode == 0, result.stderr.splitlines()[-5:]
    matches = sorted(
        glob.glob("./outputs/test_travelling_vortex/*ensemble=1*_ci_obs*.h5")
    )
    assert matches
    return matches[-1]


@pytest.mark.parametrize("da_type", ["rloc", "etpf"])
def test_da_improves_over_forecast(obs_file, da_type):
    dap = {
        "da_times": TOUT,
        "obs_attrs": OBS_ATTRS,
        "obs_path": obs_file,
        "da_type": da_type,
    }
    result = queue_run(MEMBERS, ud_payload("ci_" + da_type), dap)
    assert result.returncode == 0, result.stderr.splitlines()[-5:]

    ens_files = sorted(
        glob.glob(
            "./outputs/test_travelling_vortex/*ensemble=%i*_ci_%s*.h5"
            % (MEMBERS, da_type)
        )
    )
    assert ens_files
    t = TOUT[-1]
    with h5py.File(obs_file, "r") as truth_h5, h5py.File(ens_files[-1], "r") as ens_h5:
        for field in OBS_ATTRS:
            truth = read_field(truth_h5, field, 0, t, "after_full_step")
            rmse_fc, spread_fc = rmse_and_spread(ens_h5, truth, field, t, "before_da")
            rmse_an, spread_an = rmse_and_spread(
                ens_h5, truth, field, t, "after_full_step"
            )
            assert (
                rmse_an < rmse_fc
            ), "%s (%s): analysis RMSE %.3e !< forecast RMSE %.3e" % (
                field,
                da_type,
                rmse_an,
                rmse_fc,
            )
            assert (
                spread_an < spread_fc
            ), "%s (%s): analysis spread %.3e !< forecast spread %.3e" % (
                field,
                da_type,
                spread_an,
                spread_fc,
            )
