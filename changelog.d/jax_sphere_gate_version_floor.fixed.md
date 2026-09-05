JAX sphere full-run gates made jax-release independent. The channel and
global Williamson TC2 jax-vs-numpy checks that include the bicgstab initial
projection compared the two backends at ~1e-5, but the projection only fixes
the answer to its residual class (scipy default ``rtol=1e-5``, not overridden
by ``ud.tol``), and jax 0.11.1 lands on a different member than 0.10.1 —
~6e-4 in the channel momenta and 5.8e-4 in the global p2 with identical
numpy/scipy on the same machine — so the unpinned CI job failed on the jax
bump alone. The precision gates are now the projection-OFF
``*_stepper_bit_identical`` runs (new for the channel, hybrid + device; the
global one already existed), which agree to ~1e-16 under both releases; the
projection-ON runs remain as loose sanity checks (scalars 5e-4, momenta 2e-3,
p2 2e-3) and report all seven fields in one assertion.
