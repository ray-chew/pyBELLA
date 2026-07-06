Spherical geometry, stage B: spatially varying Coriolis —
`ud.coriolis_field` (callable of the Cartesian coordinates returning the
3 Cartesian rotation components, evaluated once per grid and cached)
feeds the H^-1 kernels and the explicit forward step as per-cell fields;
the legacy scalar `ud.coriolis_strength` path is bit-identical.
`SphericalShellMap.traditional_coriolis` provides f(phi) e_r; the JAX
backends fast-fail on rotation fields until the sphere SWE JAX stage.
