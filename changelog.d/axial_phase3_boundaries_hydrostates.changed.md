Axial-agnosticity Phase 3: hydrostatic state, vertical profiles, and
boundary handling are axis-generic — `States` profiles carry their vertical
axis (`expand_profile`), `integrated_state`/`analytical_state` integrate
along the configured vertical, the nodal-divergence wall zeroing loops over
all WALL/RAYLEIGH axes (fixing the previously broken x-WALL elliptic path),
the gravity ghost-cell handler threads the physical vertical axis (also
fixing a latent 3D bug where it read `gravity_strength[2] = 0` during
advection sweeps), `_set_boundary` mirrors the wall-normal momentum of the
actual wall axis, and the quasi-2D nodal broadcast generalises to any
degenerate axis. Bit-identical for all existing cases (verified at tol=0).
