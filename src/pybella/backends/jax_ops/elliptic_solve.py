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

import functools

import numpy as np
import jax
import jax.numpy as jnp


@functools.partial(jax.jit, static_argnames=("maxiter",))
def _solve(A, b, atol, rtol, maxiter):
    A_flat = lambda x: jnp.reshape(A(x), x.shape)
    x, _ = jax.scipy.sparse.linalg.bicgstab(
        A_flat, b, tol=rtol, atol=atol, maxiter=maxiter
    )
    return x


def bicgstab(matvec, b, atol, maxiter, rtol=1e-5):
    """Solve ``A x = b`` where ``A`` is the (already preconditioned) lap
    matvec from ``jax_ops.laplacian.lap2D/lap3D.get_linop``.

    ``matvec`` may return the boxed array (lap3D convention); it is
    re-flattened to match the solve vector. Returns a numpy float64 vector
    like the scipy path.

    get_linop returns a ``jax.tree_util.Partial`` (stable function identity,
    coefficient arrays as pytree leaves), so repeated solves at one
    resolution hit ``_solve``'s jit cache. A fresh closure per call would
    re-trace bicgstab's while_loop every step; each compiled loop pins its
    device buffers in the executable cache, which exhausts GPU memory over
    a long window (observed as CUDA_ERROR_OUT_OF_MEMORY at 1024^2).
    """
    b = jnp.asarray(b, dtype=jnp.float64)
    if isinstance(matvec, jax.tree_util.Partial):
        return np.asarray(_solve(matvec, b, atol, rtol, maxiter))
    # opaque callable: traced fresh on every call (ad-hoc operators in tests)
    A = lambda x: jnp.reshape(matvec(x), x.shape)
    x, _ = jax.scipy.sparse.linalg.bicgstab(A, b, tol=rtol, atol=atol, maxiter=maxiter)
    return np.asarray(x)
