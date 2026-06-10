Axial-agnosticity Phase 4: the 3D elliptic operator now carries the full
H^-1 tensor coefficients (C_ij = (Gamma^-1 P Theta) h[role(i),role(j)], the
same H^-1 the momentum correction applies), replacing the legacy ad-hoc x-z
`corrf` cross terms; the 3D preconditioner uses the C_ii diagonals. With
identity H^-1 the operator reduces bit-exactly to the previous one (Oracle
A), and with full Coriolis a y-uniform 3D solve now matches the trusted 2D
path under the pseudovector axis mapping to ~1e-6 (new
test_scripts/test_3d_coriolis_oracle.py). That oracle also exposed and
fixed the implicit-side sibling of the 2D out-of-plane Coriolis defect (the
pressure-correction w-row was ndim==3-only). Regenerated golden masters:
travelling_vortex_3d_coriolis (operator upgrade), igw_baldauf_brdar,
internal_long_wave, unstable_lamb (w-row fix; the igw analytic-oracle
out-of-plane error improves 0.060 -> 0.045). All other cases bit-identical.
