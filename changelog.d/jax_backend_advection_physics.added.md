JAX backend, components 3+4 of 4 (advection + physics kernels; full-run
validation): the per-sweep advection kernel (recovery + HLL upwind fluxes),
the advective rhoY flux convolution, the Coriolis H^-1 apply/coefficients,
and explicit diffusion now run on JAX under `ud.backend = "jax"`. The seam
sits inside the numpy drivers (sweeps, flips, ghost fills and flux-difference
updates are shared between backends), so flux-container semantics are
bit-faithful. New env var `PYBELLA_BACKEND=jax` flips the backend for
unmodified cases. Validation: driver-level equivalence on five states at
magnitude-scaled 1e-12 (pure explicit arithmetic), and — end to end — the
full regression suite runs on the JAX backend against the *same stored
golden-master targets* as numpy (all fields ~1e-7 max-abs vs the 1e-5
tolerance). Init-time code (hydrostatics, initial_pressure) intentionally
stays numpy: it runs once and feeds both backends identically.
