JAX backend, component 1 of 4 (pure operators): `pybella.backends.jax_ops`
mirrors `utils/operators` module-for-module (finite_difference, gradient,
convolution, divergence, preconditioner, lap2D, lap3D) with functional,
jit-clean twins of every kernel, validated against the numpy/numba
implementation as a golden master at 1e-13..1e-14 (x64, magnitude-scaled)
on realistic coefficient fields covering periodic/WALL/atmosphere
boundaries and the terrain metric. jax is an optional dependency
(`pip install pybella[jax]`); the canonical numpy path is untouched and the
equivalence suite skips cleanly without jax. New CI job `jax-equivalence`
runs the suite on Python 3.12 with an unpinned jax.
