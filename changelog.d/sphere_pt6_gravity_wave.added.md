Spherical geometry, stage D5: DCMIP-31-style nonhydrostatic gravity wave
on the small-planet compressible shell (`test_sphere_gw`) — a
latitude-independent, longitude-periodic potential-temperature
perturbation on the isothermal (Baldauf-Brdar) background whose
large-radius equatorial slice reduces to the planar B&B channel. Fixes
the general (curved-metric) free-slip gravity wall to be
FREESTREAM/WELL-BALANCED: the wall ghost now reflects the rhoY flux
`Y*(N_v.m)` as exactly odd (using the image cell's `N_v.e_up`, not the
source cell's), so the wall mass flux cancels to roundoff and J-weighted
mass/energy are conserved to machine precision under dynamics — where the
prior up-velocity reflection left an O(dz^2) leak. Confined to the
`vertical_line=False` path; vertical-line terrain and uniform-Cartesian
walls stay bit-identical (all terrain/shell/SWE golden masters green).
Stage D physics gates complete: J-weighted conservation machine-exact,
2nd-order self-convergence (Richardson 3.92), and the large-radius limit
converging to the planar Baldauf-Brdar linear oracle (zonal-field rel-L2
0.187 -> 0.098 -> 0.051 as X = 125 -> 62.5 -> 31.25).
