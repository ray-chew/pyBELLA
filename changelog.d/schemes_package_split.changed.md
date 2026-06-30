Split the 793-line `interfaces/dynamics_blending/schemes.py` god-module into a
`schemes/` package: `blending` (the `Blend` interface), `comp_psinc`, `swe_lake`,
`hydro_nonhydro` (commented hydro block kept verbatim as the reinstatement
reference) and `orchestration` (the per-timestep blending calls). Public names are
re-exported from `schemes/__init__.py`, so `time_update.py` and `prepare.py` call
sites are unchanged. Mechanical move only — FULL repro gate bit-identical at tol 0
(all 11 golden masters, incl. the blending warm bubble).
