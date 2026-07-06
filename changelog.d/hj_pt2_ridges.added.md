Hughes & Jablonowski (2023) mountain baroclinic wave, pt 2: the two
midlatitude ridges (the wave TRIGGER) via a terrain-following spherical map
plus the well-balanced "adjusted" background (`tests/test_hj_baroclinic_ridges.py`,
smoke gate `test_scripts/test_hj_baroclinic_ridges.py`). The ridges (Eq. 1,
h0 = 2000 m at 72 E / 140 E, 45 N) enter as GEOMETRY: a
`SphericalTerrainMap` (Gal-Chen radial coordinate, orography nondimensionalised
to the map's `z / h_ref` units) with no SWE bottom-topography source terms.
The balance is the surface-pressure adjustment (Eq. 2), which is exactly the
Ullrich pressure profile (Eq. B4) sampled at the surface height
`p_s(lambda, phi) = pressure(phi, z_s)`; more generally the whole adjusted
state is the analytic base state sampled at the terrain-following height
`z = r - a = Z(eta, h(lambda, phi))`. pt 1's `sol_init` is written entirely in
`metric.height`, so it reproduces the adjusted state unchanged once the map
carries the ridges — the case is pt 1 with the map swapped and the radial
axis carrying eta in [0, depth] (verified: ridge-tip surface pressure 779 hPa
vs the paper's ~773 hPa; J = r^2 cos(phi) stays positive at the ridge
latitudes only once the orography is nondimensionalised). The adjusted state
is well-balanced but not PERFECTLY so, and the residual near the ridges is
the intended trigger. Gate: over a coarse ~40 min run the ridges drive a
meridional-wind response ~8x the flat-background pt-1 adjustment (1.45 vs
0.18 m/s) LOCALISED at the ridge centres (peaks at 73 / 141 E), with the jet
bounded and nothing blowing up — the decisive ridge-triggered-initiation
signature. Next: the multi-day full-planet device run (pt 3).
