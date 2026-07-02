Added the MWR-2022 OSSE reproduction tooling: `run_scripts/osse_mwr2022.py`
drives the paper's five-run recipe per case (obs, truth, noda ensemble, LETKF
with/without blending, plus ETPF variants) with fully seeded, regenerable
observation files, and `run_scripts/osse_diagnostics.py` computes the
ensemble-mean RMSE vs truth and ensemble spread per field over the
assimilation window (CSV + PNG). The warm-bubble case gains the paper's seeded
ensemble perturbation (`delth += 10*rand()`, truth seed 1234) — inert for the
default deterministic run — and the travelling-vortex `logging.info` misuse
in `sol_init` is fixed.
