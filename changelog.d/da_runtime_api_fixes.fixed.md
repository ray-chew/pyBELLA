Fixed three run-blockers surfaced by the first N=2 assimilation cycles since
the DA repair: `dask.diagnostics` is imported explicitly (dask no longer
auto-imports it), the LETKF's scipy bindings are pinned to the intended
objects (`scipy.sparse.linalg.spsolve`, dense `scipy.linalg.eigh`,
`scipy.sparse.eye/diags` — the reference's `import scipy.sparse as sp` only
resolved under pre-1.8 scipy), and `prepare_rloc` is now built for every
`da_type` at N>1 because `obs_noiser` needs its cell/node attribute partition
(the reference had the same latent NameError for ETPF with observation noise).
LETKF-rloc and ETPF both complete a 2-cycle N=2 travelling-vortex OSSE.
