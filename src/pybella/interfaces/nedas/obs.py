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

3D conventions (there is NO native 3D pipeline, so these are the defining
conventions, gated by test_scripts/test_nedas_obs3d_selfcheck.py):

- Volumetric sparse mask, seed 778 (its own chain; the 2D seed-777
  machinery above is untouched): default_rng(778) -> per-time integer
  seeds -> per time a shuffled 0/1 mask over the inner (iicx, iicy, iicz)
  cells with ceil(N*obs_frac) zeros (observed), shared by all attrs (all
  cell-grid; p2_nodes is not observed in 3D).
- Positions under the 3D axis mapping (model.py): NEDAS x := pyBELLA
  x, NEDAS y := pyBELLA z (horizontal), NEDAS z := pyBELLA y cell-centre
  heights (vertical).
- err_std: the same VarCov reduction over the observed inner points of the
  clean 3D truth field, then max(sd, err_std_floor) — the floor guards the
  degenerate case (a z-uniform truth keeps rhow ~ 0 for all time).
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
    err_std_floor: float

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

    def _grid_of(self, attr):
        """The pyBELLA grid an attribute lives on (node for p2, else cell)."""
        return self._model.node if attr == "p2_nodes" else self._model.elem

    def _native_obs(self):
        if self._cache is not None:
            return self._cache
        if self._model.elem.ndim == 3:
            self._cache = self._native_obs_3d()
            return self._cache

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

        # sparsity masks: the exact native seed chain (1 = excluded). One
        # seed per TIME (sparse_obs_by_attr=False); each attr rebuilds the
        # mask from that seed with its OWN grid size, so cell attrs share a
        # mask and p2_nodes gets the node-grid one (utils.sparse_obs_selector).
        np.random.seed(777)
        seeds = np.random.randint(10000, size=(len(times), 1))
        if len(seeds) > 1:
            seeds = seeds.squeeze()
        masks = {}
        for tt in range(len(times)):
            for attr in self.obs_attrs:
                grid = self._grid_of(attr)
                nx, ny = grid.iicx, grid.iicy
                key = (tt, nx, ny)
                if key in masks:
                    continue
                n_obs = int(np.ceil(nx * ny * self.obs_frac))
                np.random.seed(seeds[tt])
                mask = np.array([0] * n_obs + [1] * (nx * ny - n_obs))
                np.random.shuffle(mask)
                masks[key] = mask.reshape(nx, ny)

        def mask_for(tt, attr):
            grid = self._grid_of(attr)
            return masks[(tt, grid.iicx, grid.iicy)]

        # VarCov error std per attr: per-time stds averaged over times.
        # Bit-exact native replication (utils.obs_noiser): fill a
        # (n_times, n_attrs) array and mean over axis 0 — the summation
        # order of the strided mean differs from a plain list mean by 1 ULP.
        std_dev = np.zeros((len(times), len(self.obs_attrs)))
        for tt, t in enumerate(times):
            for ai, attr in enumerate(self.obs_attrs):
                inner = clean[(t, attr)][self._grid_of(attr).i2]
                value = np.ma.array(inner, mask=mask_for(tt, attr))
                var = self.noise_percentage * ((value - value.mean()) ** 2).mean()
                std_dev[tt, ai] = var**0.5
        mean_sd = std_dev.mean(axis=0, keepdims=True)
        sd = {attr: mean_sd[0, ai] for ai, attr in enumerate(self.obs_attrs)}

        # obs positions: inner grid points (cell centres / nodes) at mask == 0
        coords = {}
        for attr in self.obs_attrs:
            grid = self._grid_of(attr)
            x1d = grid.x[grid.igx : -grid.igx]
            y1d = grid.y[grid.igy : -grid.igy]
            coords[attr] = np.meshgrid(x1d, y1d, indexing="ij")  # field-aligned

        cache = {}
        for tt, t in enumerate(times):
            time_dt = self.c.config.time_start + t * self._dt1h()
            for attr in self.obs_attrs:
                sel = mask_for(tt, attr) == 0
                xg, yg = coords[attr]
                values = clean[(t, attr)][self._grid_of(attr).i2][sel]
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

    def _native_obs_3d(self):
        """Volumetric 3D obs (conventions in the module docstring)."""
        times = [round(float(t), 3) for t in self.da_times]
        elem = self._model.elem
        nx, ny, nz = elem.iicx, elem.iicy, elem.iicz

        clean = {}
        with h5py.File(self.obs_file, "r") as f:
            for t in times:
                for attr in self.obs_attrs:
                    label = "%s_ensemble_mem=0_%.3f_after_full_step" % (attr, t)
                    clean[(t, attr)] = np.squeeze(f[attr][label][:])

        # seed-778 volumetric masks: one per time, shared by all attrs
        rng = np.random.default_rng(778)
        time_seeds = rng.integers(0, 2**31 - 1, size=len(times))
        n_pts = nx * ny * nz
        n_obs = int(np.ceil(n_pts * self.obs_frac))
        masks = []
        for tt in range(len(times)):
            rng_t = np.random.default_rng(int(time_seeds[tt]))
            mask = np.array([0] * n_obs + [1] * (n_pts - n_obs))
            rng_t.shuffle(mask)
            masks.append(mask.reshape(nx, ny, nz))

        # VarCov err std over the OBSERVED inner points, then the floor
        std_dev = np.zeros((len(times), len(self.obs_attrs)))
        for tt, t in enumerate(times):
            for ai, attr in enumerate(self.obs_attrs):
                inner = clean[(t, attr)][elem.i2]
                value = np.ma.array(inner, mask=masks[tt])
                var = self.noise_percentage * ((value - value.mean()) ** 2).mean()
                std_dev[tt, ai] = var**0.5
        mean_sd = std_dev.mean(axis=0, keepdims=True)
        floor = float(getattr(self, "err_std_floor", 0.0) or 0.0)
        sd = {
            attr: max(mean_sd[0, ai], floor) for ai, attr in enumerate(self.obs_attrs)
        }

        # positions: inner cell centres under the 3D axis mapping (model.py)
        x1d = elem.x[elem.igx : -elem.igx]
        y1d = elem.y[elem.igy : -elem.igy]
        z1d = elem.z[elem.igz : -elem.igz]
        xg, yg, zg = np.meshgrid(x1d, y1d, z1d, indexing="ij")  # field-aligned

        cache = {}
        for tt, t in enumerate(times):
            time_dt = self.c.config.time_start + t * self._dt1h()
            sel = masks[tt] == 0
            for attr in self.obs_attrs:
                values = clean[(t, attr)][elem.i2][sel]
                nobs = values.size
                cache[(t, attr)] = {
                    "obs": values,
                    "t": np.full(nobs, time_dt),
                    "x": xg[sel],  # pyBELLA x  -> NEDAS x (horizontal)
                    "y": zg[sel],  # pyBELLA z  -> NEDAS y (horizontal)
                    "z": yg[sel],  # pyBELLA y  -> NEDAS z (vertical, centres)
                    "err_std": np.full(nobs, sd[attr]),
                }
        return cache

    @staticmethod
    def _dt1h():
        from NEDAS.utils.conversion import dt1h

        return dt1h

    def read_obs_from_file(self, **kwargs):
        raise NotImplementedError("read_obs serves the native obs file directly")
