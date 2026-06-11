"""JAX elliptic solve: bicgstab with scipy-matching convergence semantics.

Component 2 of the JAX migration. The seam sits at the linear-solver level
(`implicit_euler.do_implicit_part`): host-side coefficient assembly (ghost
cells, operator coefficients, terrain folding, diagonal preconditioning,
rhs construction) stays on the canonical numpy path, so both backends solve
the same preconditioned system (bit-identical up to the ulp-level Coriolis
coefficients, which component 3 routed through the JAX twin); the operator
application and the Krylov iteration run in JAX.

Convergence matching: scipy's ``bicgstab(A, b, atol=ud.tol)`` stops at
``||r|| <= max(rtol * ||b||, atol)`` with ``rtol = 1e-5`` (scipy default);
``jax.scipy.sparse.linalg.bicgstab(tol=rtol, atol=atol)`` implements the
same criterion. The iterates themselves differ between implementations, so
solutions agree to solver tolerance, not bitwise — equivalence is asserted
via the residual certificate (the JAX solution satisfies scipy's stopping
criterion measured with the *numpy* operator) plus field-level agreement.
"""

import numpy as np
import jax
import jax.numpy as jnp


def bicgstab(matvec, b, atol, maxiter, rtol=1e-5):
    """Solve ``A x = b`` where ``A`` is the (already preconditioned) lap
    matvec from ``jax_ops.laplacian.lap2D/lap3D.get_linop``.

    ``matvec`` may return the boxed array (lap3D convention); it is
    re-flattened to match the solve vector. Returns a numpy float64 vector
    like the scipy path.
    """
    b = jnp.asarray(b, dtype=jnp.float64)
    A = lambda x: jnp.reshape(matvec(x), x.shape)
    x, _ = jax.scipy.sparse.linalg.bicgstab(A, b, tol=rtol, atol=atol, maxiter=maxiter)
    return np.asarray(x)
