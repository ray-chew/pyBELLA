General curvilinear gradient map and elliptic fold (tfc pt 3):
`apply_gradient_map` now applies A_{ka} = (N_a)_k / J and
`elliptic_tensor(_2d)` folds M = (1/J) N H^-1 N^T (new
`elliptic_diag_geometric` for the 2D preconditioner diagonal, H^-1 still
excluded). Bit-exact h == 0 reduction, one-ulp vertical-line reduction;
laplacian kernels untouched; the JAX device path shares the same functions.
Oracle suites extended with genuinely stretched Tier-2 maps (3D + native
2D + Coriolis) and an SPD/symmetry gate.
