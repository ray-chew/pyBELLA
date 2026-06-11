JAX backend, device-resident plan phase A: functional, jitted twins of every
per-step boundary operation — ghost-cell fills (periodic/WALL no-gravity
pads, the sequential two-layer hydrostatic gravity fill incl. terrain
contravariant reflection and ATMOSPHERIC_EXTENSION), ghost-node fills (the
periodic overlap exchange, reflect, quasi-2D degenerate broadcast), Rayleigh
sponge damping (with func-mode forcing arrays), and wall-node rhs scaling —
routed through the hybrid seams under `ud.backend = "jax"`. All gravity-fill
index algebra, stratification-at-ghost-coordinates and metric slices are
precomputed per orientation (canonical and sweep-flipped) at config build;
the kernels are pure gathers + closed-form arithmetic. Validated per-BC on
six states × flipped sweep orientations × sol-override ×
compressible/incompressible at magnitude-scaled 1e-13, with np.pad-callback
micro-tests (bit-faithful including the degenerate quasi-2D axis), plus the
full 10-case regression gate on the hybrid backend. After this phase the
hybrid step has no numpy kernels left — only orchestration.
