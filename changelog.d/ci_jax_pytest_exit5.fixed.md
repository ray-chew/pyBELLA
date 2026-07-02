CI: drop the redundant `test_jax_*` pytest invocations from the numpy
integration job — without jax they collect zero tests, which pytest >= 9
fails with exit code 5; the jax-equivalence job runs them.
