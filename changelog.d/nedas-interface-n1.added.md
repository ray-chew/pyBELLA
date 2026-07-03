Phase N1 of the NEDAS interface: working `PyBellaModel`/`PyBellaObs` adapter
(in-memory per-member ModelState forecasts, native-parity observations),
`run_scripts/nedas_run.py` launcher + `nedas_osse_diagnostics.py` exporter,
and `test_scripts/test_nedas_obs_parity.py` (obs byte-identical to the frozen
native pipeline). TV noda/EnDA/EnDAB OSSEs run end-to-end through NEDAS.
