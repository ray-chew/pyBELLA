Spherical geometry, JAX twins (task 1a): the hybrid JAX backend now fills
ghost cells for non-vertical-line (spherical) metrics.
`backends/jax_ops/boundary.py` gains (1) the general free-slip WALL mirror
(twin of `cell_boundary._mirror_momenta_general`: symmetric pad + per-side
contravariant Cramer reflection with the local area normals), invoked at
canonical orientation for the phi and thin-shell degenerate-r walls, and
(2) the well-balanced radial gravity ghost fill (twin of the e_up branch of
`_calculate_ghost_values`/`_assign_ghost_values`: hydrostatic rho/rhoY,
tangential-velocity copy, and the beta·e_up momentum reconstruction that
makes the rhoY wall flux exactly odd), in phys + sweep orientations. The
non-vertical-line fast-fail guards are removed. Gated by
`test_cells_sphere_*` in `test_jax_boundary_equivalence.py` and by
`test_strang_sphere_{swe,gw}` (full sweep driver) in
`test_jax_advection_equivalence.py` — JAX vs numpy at the ulp floor,
including the mid-sweep (flipped) orientations.
