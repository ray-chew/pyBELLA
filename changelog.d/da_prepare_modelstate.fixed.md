Rebuilt `data_assimilation/prepare.py` on the ModelState/EnsembleState API:
ensemble members are constructed from fresh `CellSolField`/`NodePressureField`
containers with one seeded `sol_init` call each (matching the paper-era
reference, where `sol_init` never pre-ran for N>1 — re-running it on the
initialised member 0 double-applies the `+=` initialisations and blows up the
forecast), each member gets its own `FlowSolverCache`, and ghost cells are set
at construction. Also: observation loading is skipped when `da_times` is empty
(pure ensemble forecasts need no obs file), `obs_path` now defaults to `None`
with an actionable error pointing at the dap-rewrite route, and the broken
`es.flux`/`ensembble_state` references are gone. Verified with an N=2
travelling-vortex ensemble forecast smoke run.
