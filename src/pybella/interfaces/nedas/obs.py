"""PyBellaObs — serve the native OSSE observations to NEDAS, byte-identically.

Subclasses ``NEDAS.core.Dataset`` directly, NOT SyntheticObs: any SyntheticObs
instance makes the scheme generate its own truth and add unseeded noise
(NEDAS/core/obs.py collect_obs_seq). The read_obs branch instead serves the
native observation HDF5 (produced by ``run_scripts/osse_mwr2022.py --runs
obs``) verbatim.

Native conventions replicated (data_assimilation is frozen and NOT imported;
verified against src/pybella/data_assimilation/{params,utils,analysis}.py):

- The analyses assimilate the CLEAN obs fields — obs_noisy is computed but
  never used; only the noise model enters, through the obs error covariance.
  So obs values here are the H5 fields at the sparse points, no noise added.
- Sparsity mask: seed 777 -> randint(10000, size=(n_times, 1)).squeeze(),
  one seed per analysis time (sparse_obs_by_attr=False: all attrs share the
  time's mask); per time, K=ceil(N*obs_frac) zeros (observed) shuffled into
  an (iicx, iicy) array. (utils.sparse_obs_selector)
- Error std (VarCov): var = noise_percentage * mean((v - mean(v))^2) over the
  UNMASKED inner points of the clean field, per attr per time, then the std
  is AVERAGED over times per attr. (utils.obs_noiser)
- H5 labels: /<attr>/<attr>_ensemble_mem=0_<t:.3f>_after_full_step, full
  arrays including ghost frames. (params.init.load_obs)
"""

import h5py
import numpy as np
from NEDAS.core import Dataset


class PyBellaObs(Dataset):
    """Native pyBELLA OSSE observations as a NEDAS (non-synthetic) dataset.

    Config keys (dataset_def.pybella):
        obs_file: path to the native obs HDF5
        obs_attrs: attributes stored in the obs file, native order
        da_times: the native analysis times (defines the seed-777 mask chain)
        obs_frac: fraction of inner points observed (native default 0.10)
        noise_percentage: VarCov fraction (native default 0.05)
    """

    obs_file: str
    obs_attrs: list
    da_times: list
    obs_frac: float
    noise_percentage: float

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if "model_src" in kwargs:
            self._model = self.c.models[kwargs["model_src"]]
            self.variables = dict(self._model.variables)
            self.grid = self._model.grid
        self._cache = None

    def read_obs(self, **kwargs):
        """Return the obs_seq for (cycle time, variable) from the native file."""
        kwargs = self.parse_kwargs(kwargs)
        t = self._model.nondim_time(kwargs["time"])
        seq = self._native_obs()[(t, kwargs["name"])]
        return {
            key: value.copy() if isinstance(value, np.ndarray) else value
            for key, value in seq.items()
        }

    # --- native pipeline replication -------------------------------------

    def _native_obs(self):
        if self._cache is not None:
            return self._cache

        elem = self._model.elem
        i2 = elem.i2
        nx, ny = elem.iicx, elem.iicy
        npts = nx * ny
        times = [round(float(t), 3) for t in self.da_times]

        # clean fields from the native obs H5
        clean = {}
        with h5py.File(self.obs_file, "r") as f:
            for t in times:
                for attr in self.obs_attrs:
                    label = "%s_ensemble_mem=0_%.3f_after_full_step" % (attr, t)
                    data = f[attr][label][:]
                    if data.ndim == 3:  # 3D run: horizontal slice, as native
                        data = data[:, 0, :]
                    clean[(t, attr)] = data

        # sparsity masks: the exact native seed chain (1 = excluded)
        np.random.seed(777)
        seeds = np.random.randint(10000, size=(len(times), 1))
        if len(seeds) > 1:
            seeds = seeds.squeeze()
        n_obs = int(np.ceil(npts * self.obs_frac))
        masks = {}
        for tt in range(len(times)):
            np.random.seed(seeds[tt])
            mask = np.array([0] * n_obs + [1] * (npts - n_obs))
            np.random.shuffle(mask)
            masks[tt] = mask.reshape(nx, ny)

        # VarCov error std per attr: per-time stds averaged over times.
        # Bit-exact native replication (utils.obs_noiser): fill a
        # (n_times, n_attrs) array and mean over axis 0 — the summation
        # order of the strided mean differs from a plain list mean by 1 ULP.
        std_dev = np.zeros((len(times), len(self.obs_attrs)))
        for tt, t in enumerate(times):
            for ai, attr in enumerate(self.obs_attrs):
                value = np.ma.array(clean[(t, attr)][i2], mask=masks[tt])
                var = self.noise_percentage * ((value - value.mean()) ** 2).mean()
                std_dev[tt, ai] = var**0.5
        mean_sd = std_dev.mean(axis=0, keepdims=True)
        sd = {attr: mean_sd[0, ai] for ai, attr in enumerate(self.obs_attrs)}

        # obs positions: inner cell centres where mask == 0
        x1d = elem.x[elem.igx : -elem.igx]
        y1d = elem.y[elem.igy : -elem.igy]
        xg, yg = np.meshgrid(x1d, y1d, indexing="ij")  # (nx, ny), field-aligned

        cache = {}
        for tt, t in enumerate(times):
            sel = masks[tt] == 0
            time_dt = self.c.config.time_start + t * self._dt1h()
            for attr in self.obs_attrs:
                values = clean[(t, attr)][i2][sel]
                nobs = values.size
                cache[(t, attr)] = {
                    "obs": values,
                    "t": np.full(nobs, time_dt),
                    "z": np.zeros(nobs),
                    "y": yg[sel],
                    "x": xg[sel],
                    "err_std": np.full(nobs, sd[attr]),
                }
        self._cache = cache
        return cache

    @staticmethod
    def _dt1h():
        from NEDAS.utils.conversion import dt1h

        return dt1h

    def read_obs_from_file(self, **kwargs):
        raise NotImplementedError("read_obs serves the native obs file directly")
