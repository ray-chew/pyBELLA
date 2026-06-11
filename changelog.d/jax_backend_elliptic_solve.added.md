JAX backend, component 2 of 4 (elliptic solve): new `ud.backend` flag
(default `"numpy"`) selects the solver for the semi-implicit pressure
projection. With `backend = "jax"`, `implicit_euler.do_implicit_part` builds
the JAX laplacian linop and solves with `jax.scipy.sparse.linalg.bicgstab`
under scipy-matching convergence semantics (rtol 1e-5 / atol `ud.tol`);
coefficient assembly, preconditioning and the rhs stay on the numpy path, so
both backends solve bit-identical systems. Validated by residual
certificates (the JAX solution satisfies scipy's stopping criterion under
the numpy operator) and end-to-end `do_implicit_part` equivalence on six
states at magnitude-scaled 1e-5 (the convergence-slack floor; the vortex
case agrees at 3.6e-10). Default-path behaviour is unchanged.
