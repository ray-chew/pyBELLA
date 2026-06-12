General curvilinear divergence fluxes (tfc pt 2): the terrain rhs now
assembles `F_a = N_a . (theta m)` (vertical-first contraction, bit-exact
reduction to the legacy J-weighted/contravariant fluxes) in both the numpy
and JAX backends; differencing stencils untouched. New freestream-preservation
tripwire `test_scripts/test_freestream.py` (uniform flow on a wavy 3D Tier-2
map: defect bounded and second-order convergent; bit-exact zero on flat).
