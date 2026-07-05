NEDAS Phase E: the first true 3D DA — igw3d OSSE through NEDAS. The adapter
gains a 3D branch (horizontal (x,z) grid + vertical levels with per-level
read/write and interface z_coords), PyBellaObs a seed-778 volumetric obs
convention with a VarCov err-std floor, the internal-long-wave IC a 3D
branch + seeded theta'-wave member perturbations (2D path bit-identical;
oracle-gated), and the OSSE tooling igw3d generation, ndim-general
diagnostics with the w=rhov/rho blending probe, and dev/prod configs.
Self-anchored gates (no native 3D pipeline exists):
test_igw3d_ic_oracle.py, test_nedas_obs3d_selfcheck.py, and the igw3d
variant of test_nedas_vmap_gate.py. Findings + ladder results in
dev_notes/nedas_interface.md Phase E.
