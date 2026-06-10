Terrain phase 6: advection uses contravariant vertical / J-weighted horizontal
mass fluxes (HLL upwinding follows automatically), the cell update divides by
J, recovery's Courant velocity divides J back out, and the CFL gains the
metric vertical signal speed. Orography is evaluated on periodic-wrapped
coordinates — a non-periodic hill on a periodic axis otherwise makes the
elliptic system inconsistent at the duplicated nodes (found via dense
null-space analysis; bicgstab diverged). smoke_agnesi now runs the 400 m
witch-of-Agnesi hill end-to-end: J-weighted mass/P conserved to 1e-10,
max |w| ≈ 0.44 m/s vs the ~0.5 m/s linear estimate.
