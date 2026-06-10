Added a physics oracle for the Baldauf & Brdar (2013) internal-gravity-wave
case: `pybella/tests/baldauf_brdar_analytic.py` builds a numerically-exact
solution of the linearised compressible Euler equations about the isothermal
background (Bretherton-transformed constant-coefficient system, staggered
vertical collocation, exact-in-time eigenpropagation per x-Fourier mode;
energy drift ~1e-13) and evolves the simulation's own initial condition.
`test_scripts/test_igw_analytic.py` gates the regression configuration
against it (catching sign/dispersion/amplitude *wrongness*, complementing
the golden masters which catch *change*) and writes ref/sim/diff PNGs.
Refinement studies (dt 500->125 s, dx 20->10 km, f on/off) decompose the
measured sim-vs-linear residual and surface a known solver limitation: in 2D
runs the out-of-plane momentum receives only the implicit half of the
Coriolis rotation (`explicit_euler.do_forward_step` skips the `rhow` row for
`ndim == 2`), pinning that component's error at ~0.44 rel-L2.
