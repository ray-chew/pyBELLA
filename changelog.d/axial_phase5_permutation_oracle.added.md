Axial-agnosticity Phase 5: the permutation oracle
(`test_scripts/test_permutation_oracle.py`, CI-wired) proves the endgame —
the 2D internal-long-wave reference (gravity, stratification, walls, full
Coriolis) embedded as z-vertical (gravity_direction=2) and x-vertical
(gravity_direction=0) quasi-2D 3D twins reproduces the sigma-mapped 2D
fields through full solver steps to <=1e-6, with exact uniformity along
the degenerate axis and the Coriolis pseudovector mapping produced
automatically by the role-based configuration. Also fixes a latent bug the
oracle exposed: SpaceDiscr stored `dxyz/ig/ic/stride` as class-level shared
arrays, so two grids coexisting in one process corrupted each other; they
are now per-instance. Bit-identical for all existing cases.
