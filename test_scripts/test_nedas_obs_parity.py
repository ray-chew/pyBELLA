"""Byte-parity check: PyBellaObs vs the frozen native obs pipeline.

PyBellaObs REIMPLEMENTS the native sparse-mask (seed 777) and VarCov err-std
conventions rather than importing the frozen ``data_assimilation`` layer; this
test is the only place the two implementations meet. It asserts, per analysis
time and attribute, that

* the observed-point set (mask == 0) is identical,
* the obs values equal the clean H5 fields at those points (both frameworks
  assimilate the CLEAN obs; noise enters only via the error covariance),
* PyBellaObs's err_std matches sqrt(native obs_covar) exactly.

Needs the regenerable TV obs file (run_scripts/osse_mwr2022.py tv --runs obs);
exits 0 with a SKIP message when it is missing.

Run: python test_scripts/test_nedas_obs_parity.py
"""

import glob
import os
import sys

import numpy as np

OBS_GLOB = "outputs/test_travelling_vortex/*ensemble=1_64_64*_obs.h5"
DA_TIMES = [round(float(t), 3) for t in np.arange(0.25, 3.25, 0.25)]
OBS_ATTRS = ["rhou", "rhov"]
K = 10


def native_pipeline(obs_file):
    """Run the frozen data_assimilation obs prep (test-only import)."""
    import importlib

    from pybella.data_assimilation import letkf as da_letkf
    from pybella.data_assimilation import params as da_params
    from pybella.data_assimilation import utils as da_utils
    from pybella.flow_solver.discretisation import grid as dis_grid
    from pybella.interfaces.ic_config import IC_MODULES
    from pybella.utils import user_data

    case = importlib.import_module(IC_MODULES["test_travelling_vortex"])
    ud = user_data.UserDataInit(**vars(case.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)

    dap = da_params.init(K, da_type="rloc")
    dap.update_dap({"da_times": DA_TIMES, "obs_attrs": OBS_ATTRS, "obs_path": obs_file})
    obs = dap.load_obs(dap.obs_path)
    rloc = da_letkf.prepare_rloc(ud, elem, node, dap, K)
    mask = da_utils.sparse_obs_selector(obs, elem, node, ud, dap)
    _, obs_covar = da_utils.obs_noiser(obs, mask, dap, rloc, elem)
    return elem, obs, mask, obs_covar, rloc


def nedas_obs_seqs(obs_file):
    """Instantiate the NEDAS context from the run config; pull all obs_seqs."""
    import pybella.interfaces.nedas as pybella_nedas

    pybella_nedas.register()
    from NEDAS.core.context import Context
    from NEDAS.utils.conversion import dt1h

    c = Context(config_file="run_scripts/nedas_tv_osse.yml", nens=K)
    dataset = c.datasets["pybella"]
    dataset.obs_file = obs_file  # in case the config points elsewhere
    seqs = {}
    for t in DA_TIMES:
        for attr in OBS_ATTRS:
            seqs[(t, attr)] = dataset.read_obs(
                name=attr, time=c.config.time_start + t * dt1h
            )
    return seqs


def main():
    matches = sorted(glob.glob(OBS_GLOB))
    if not matches:
        print("SKIP: no obs file (%s); run osse_mwr2022.py tv --runs obs" % OBS_GLOB)
        return 0
    obs_file = matches[-1]
    print("obs file:", obs_file)

    elem, obs, mask, obs_covar, rloc = native_pipeline(obs_file)
    seqs = nedas_obs_seqs(obs_file)

    i2 = elem.i2
    x1d = elem.x[elem.igx : -elem.igx]
    y1d = elem.y[elem.igy : -elem.igy]
    xg, yg = np.meshgrid(x1d, y1d, indexing="ij")

    checked = 0
    for tt, t in enumerate(DA_TIMES):
        for attr in OBS_ATTRS:
            seq = seqs[(t, attr)]
            sel = np.asarray(mask[tt][attr])[i2] == 0
            assert seq["obs"].size == sel.sum(), (t, attr, "nobs mismatch")
            assert np.array_equal(seq["obs"], np.asarray(obs[tt][attr])[i2][sel])
            assert np.array_equal(seq["x"], xg[sel]) and np.array_equal(
                seq["y"], yg[sel]
            )
            cidx = rloc.ca.index(attr)
            native_sd = float(np.sqrt(obs_covar[0][tt, cidx]))
            assert np.allclose(seq["err_std"], native_sd, rtol=0, atol=0), (
                t,
                attr,
                seq["err_std"][0],
                native_sd,
            )
            checked += 1
    print(
        "PASS: %d (time, attr) obs records byte-identical to the native pipeline"
        % checked
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
