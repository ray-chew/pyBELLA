Terrain phase 1: the nodal divergence computes `J∇·F` when terrain is active —
J-weighted horizontal fluxes and the contravariant vertical flux
`θ(mom_v − G1·mom_h1 − G2·mom_h2)` in role space, with unchanged differencing
stencils. Ghost-slab wall zeroing covers the contravariant fluxes automatically.
Forced-flat oracle (≤1e-13 vs plain path) and role-wiring checks added.
