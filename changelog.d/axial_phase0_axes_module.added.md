Added the axis-geometry module `pybella/utils/axes.py` (role-space convention:
cyclic permutation mapping (h1, v, h2) roles onto array axes, vertical-axis
accessors, slab/profile/permutation helpers) with unit tests, plus the
bit-for-bit H5 run comparator `test_scripts/compare_h5_runs.py` used to gate
the pure-refactor phases of the axial-agnosticity work. No behaviour change.
