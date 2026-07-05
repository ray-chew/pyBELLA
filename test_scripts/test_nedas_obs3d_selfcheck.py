"""Self-anchored 3D obs gate (Phase E): PyBellaObs 3D vs an independent
re-derivation of its own documented conventions.

There is NO native 3D pipeline to be byte-identical to (the frozen
data_assimilation layer is 2D-plane-only), so this gate re-derives the
seed-778 volumetric mask chain and the VarCov+floor err-std formula
in-test — without importing PyBellaObs internals — and asserts, per
(analysis time, attribute):

* the observed-point set is identical and n_obs == ceil(N * obs_frac),
* obs values are byte-equal to the clean H5 fields at those points,
* positions equal the inner cell-centre meshgrid under the documented
  axis mapping (NEDAS x/y = pyBELLA x/z horizontal, NEDAS z = pyBELLA y),
* err_std equals the re-derived VarCov std with the err_std_floor applied,
* two PyBellaObs instantiations agree (determinism).

Needs the regenerable igw3d obs file (run_scripts/osse_mwr2022.py igw3d
--runs obs); exits 0 with a SKIP message when missing.

Run: python test_scripts/test_nedas_obs3d_selfcheck.py
"""

import glob
import sys

import h5py
import numpy as np

OBS_GLOB = "outputs/test_internal_long_wave/*ensemble=1_65_16_16*_obs.h5"
CONFIG = "run_scripts/nedas_igw3d_enda.yml"
OBS_ATTRS = ["rho", "rhou", "rhov", "rhow", "rhoY"]
DA_TIMES = [round(float(t), 3) for t in np.arange(60.0, 660.0, 60.0)]
OBS_FRAC = 0.10
NOISE_PCT = 0.05
ERR_FLOOR = 1.0e-6
K = 10


def nedas_obs(obs_file):
    """Instantiate the NEDAS context from the run config (as the 2D parity
    test does); return (model, obs-cache dict)."""
    import pybella.interfaces.nedas as pybella_nedas

    pybella_nedas.register()
    from NEDAS.core.context import Context

    c = Context(config_file=CONFIG, nens=K)
    dataset = c.datasets["pybella"]
    dataset.obs_file = obs_file
    dataset._cache = None
    return c.models["pybella"], dataset._native_obs()


def main() -> int:
    files = sorted(glob.glob(OBS_GLOB))
    if not files:
        print("SKIP: no igw3d obs file (run osse_mwr2022.py igw3d --runs obs)")
        return 0
    obs_file = files[-1]

    model, cache = nedas_obs(obs_file)
    elem = model.elem
    nx, ny, nz = elem.iicx, elem.iicy, elem.iicz
    cache2 = nedas_obs(obs_file)[1]

    # independent re-derivation of the documented conventions
    rng = np.random.default_rng(778)
    time_seeds = rng.integers(0, 2**31 - 1, size=len(DA_TIMES))
    n_pts = nx * ny * nz
    n_obs = int(np.ceil(n_pts * OBS_FRAC))

    x1d = elem.x[elem.igx : -elem.igx]
    y1d = elem.y[elem.igy : -elem.igy]
    z1d = elem.z[elem.igz : -elem.igz]
    xg, yg, zg = np.meshgrid(x1d, y1d, z1d, indexing="ij")

    checked = 0
    with h5py.File(obs_file, "r") as f:
        std_dev = np.zeros((len(DA_TIMES), len(OBS_ATTRS)))
        sels = []
        cleans = {}
        for tt, t in enumerate(DA_TIMES):
            rng_t = np.random.default_rng(int(time_seeds[tt]))
            mask = np.array([0] * n_obs + [1] * (n_pts - n_obs))
            rng_t.shuffle(mask)
            mask = mask.reshape(nx, ny, nz)
            sels.append(mask == 0)
            for ai, attr in enumerate(OBS_ATTRS):
                label = "%s_ensemble_mem=0_%.3f_after_full_step" % (attr, t)
                clean = np.squeeze(f[attr][label][:])[elem.i2]
                cleans[(t, attr)] = clean
                value = np.ma.array(clean, mask=mask)
                var = NOISE_PCT * ((value - value.mean()) ** 2).mean()
                std_dev[tt, ai] = var**0.5
        mean_sd = std_dev.mean(axis=0, keepdims=True)

        for tt, t in enumerate(DA_TIMES):
            sel = sels[tt]
            for ai, attr in enumerate(OBS_ATTRS):
                seq = cache[(t, attr)]
                assert seq["obs"].size == n_obs, (t, attr, "n_obs mismatch")
                assert np.array_equal(seq["obs"], cleans[(t, attr)][sel]), (
                    t, attr, "obs values differ from clean field at mask points"
                )
                assert np.array_equal(seq["x"], xg[sel]), (t, attr, "x")
                assert np.array_equal(seq["y"], zg[sel]), (t, attr, "y (pyB z)")
                assert np.array_equal(seq["z"], yg[sel]), (t, attr, "z (pyB y)")
                want_sd = max(mean_sd[0, ai], ERR_FLOOR)
                assert np.all(seq["err_std"] == want_sd), (
                    t, attr, seq["err_std"][0], want_sd
                )
                s2 = cache2[(t, attr)]
                for k in ("obs", "x", "y", "z", "err_std"):
                    assert np.array_equal(seq[k], s2[k]), (t, attr, k, "nondet")
                checked += 1

    print(
        "PASS: %d (time, attr) 3D obs records match the re-derived seed-778 "
        "mask, clean values, axis-mapped positions and VarCov+floor err_std"
        % checked
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
