Hughes & Jablonowski (2023) mountain baroclinic wave, pt 1: the well-balanced
Ullrich base state on the deep spherical shell as a pyBELLA initial condition
(no topography) with a steadiness gate (`tests/test_hj_baroclinic.py`,
`test_scripts/test_hj_baroclinic.py`). The 3D balance is split the way the
discretisation carries it: the VERTICAL hydrostatic balance goes into a
z-only equatorial-Ullrich field-mode HydroState reference (built by the
terrain-quadrature branch of `hydrostatics.integrated_state` from a custom
`ud.stratification`), while the MERIDIONAL pressure structure goes into the
Exner perturbation `p2_nodes` — the momentum pressure-gradient force uses only
grad(p2), so the meridional gradient the Coriolis force on the jet must
balance has to live there (the 3D analogue of TC2's geostrophic depth in p2).
Full Coriolis is the constant embedded rotation vector `2*Omega_nd*(0,0,+1)`
(the mirrored-embedding pseudovector flip). The registered case
`test_hj_baroclinic` runs `SphericalShellMap` (full-size a, deep ~30 km shell,
lambda-periodic, radial gravity, +-80 deg latitude free-slip walls). Gate:
with no ridges the background stays steady — the ~28 m/s midlatitude jet holds
to <0.05%, the geostrophic adjustment saturates the meridional wind at
~0.18 m/s (0.7% of the jet) and then plateaus, nothing blows up; a wrong
Coriolis sign/factor would drive O(jet) meridional wind within a few steps.
Next: the ridges via `SphericalTerrainMap` + the Eq. 2 adjusted balance, then
the multi-day full-planet device run.
