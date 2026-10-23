Strict steady linear mountain-wave oracle (`tests/linear_mountain_oracle.py`)
for arbitrary periodic terrain on the model's own domain, linearised about the
model's own constant-N (non-isothermal) hydrostatic background: Boussinesq or
fully compressible, hydrostatic or nonhydrostatic, with a Liouville-Green
radiation condition. Self-tested in `test_scripts/test_linear_oracle.py`
(isothermal exactness, Boussinesq and Smith limits, Eliassen-Palm flux
constancy, PDE residual, step/radiation-height convergence). Study driver and
plots in `run_scripts/strict_linear_oracle/`.
