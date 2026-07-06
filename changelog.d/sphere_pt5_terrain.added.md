Spherical geometry, stage E: terrain-following coordinates on the sphere
— `SphericalTerrainMap` composes the existing `VerticalTransform`
(Gal-Chen) into the radial coordinate, r = a + Z(eta, h(lambda, phi)),
with lambda-wrapped orography, radial coordinate lines (t_eta || e_r)
and the new `up_direction` map hook keeping GRAVITY radial while the
coordinate surfaces tilt with the slope. Gates: h == 0 reduces to the
shell map; duality + 2nd-order metric identity over a wavy hill (the
freestream tripwire); the resting isothermal atmosphere over a spherical
Gaussian-belt hill stays at rest to the solve floor.
