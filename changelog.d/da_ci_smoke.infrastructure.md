Added `test_scripts/test_da_smoke.py` to CI: a seeded 2-cycle, N=4
travelling-vortex OSSE for both LETKF (rloc) and ETPF, asserting that the
analysis ensemble mean beats the forecast against the regenerated truth on
the observed momentum fields and that the analysis spread contracts.
Statistical, not bitwise (the dask-chunked rloc analysis makes bitwise
comparisons fragile). The DA layer is frozen at the MWR-2022 reproduction:
2D x-y, vertical=1, numpy backend, bug fixes only — superseded by the NEDAS
interface as the maintained engine.
