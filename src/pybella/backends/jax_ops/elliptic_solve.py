"""JAX elliptic solve: bicgstab with scipy-matching convergence semantics.

The seam sits at the linear-solver level (`implicit_euler.do_implicit_part`):
host-side coefficient assembly (ghost cells, operator coefficients, terrain
folding, diagonal preconditioning, rhs construction) stays on the canonical
numpy path, so both backends solve the same preconditioned system
(bit-identical up to the ulp-level Coriolis coefficients, which the JAX
Coriolis twin supplies); the operator application and the Krylov iteration
run in JAX.

Convergence matching: scipy's ``bicgstab(A, b, atol=ud.tol)`` stops at
``||r|| <= max(rtol * ||b||, atol)`` with ``rtol = 1e-5`` (scipy default);
``jax.scipy.sparse.linalg.bicgstab(tol=rtol, atol=atol)`` implements the
same criterion. The iterates themselves differ between implementations, so
solutions agree to solver tolerance, not bitwise — equivalence is asserted
via the residual certificate (the JAX solution satisfies scipy's stopping
criterion measured with the *numpy* operator) plus field-level agreement.
"""

import functools
import logging

import numpy as np
import jax
import jax.numpy as jnp


@functools.partial(jax.jit, static_argnames=("maxiter",))
def _solve(A, b, x0, atol, rtol, maxiter):
    A_flat = lambda x: jnp.reshape(A(x), x.shape)
    x, _ = jax.scipy.sparse.linalg.bicgstab(
        A_flat, b, x0=x0, tol=rtol, atol=atol, maxiter=maxiter
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
    A = lambda x: jnp.reshape(matvec(x), x.shape)
    # Restart-on-stall, mirroring implicit_euler._bicgstab_with_restarts: jax's
    # bicgstab reports no breakdown, so the achieved residual decides. The
    # sphere initial projection needs this (see the numpy twin's docstring).
    x0 = None
    for attempt in range(_MAX_RESTARTS + 1):
        if isinstance(matvec, jax.tree_util.Partial):
            x = _solve(matvec, b, x0, atol, rtol, maxiter)
        else:
            # opaque callable: traced fresh on every call (ad-hoc operators in tests)
            x, _ = jax.scipy.sparse.linalg.bicgstab(
                A, b, x0=x0, tol=rtol, atol=atol, maxiter=maxiter
            )
        if attempt == _MAX_RESTARTS or _residual_ok(A, b, x, atol, rtol):
            break
        x0 = x
    _check_residual(A, b, x, atol, rtol, maxiter)
    return np.asarray(x)


_MAX_RESTARTS = 5


def _residual_ok(A, b, x, atol, rtol):
    nb = float(jnp.linalg.norm(b))
    nr = float(jnp.linalg.norm(b - A(x)))
    return nr <= 10.0 * max(rtol * nb, atol)  # same slack as the numpy check


def _check_residual(A, b, x, atol, rtol, maxiter):
    """jax's bicgstab returns no ``info``: a breakdown or a maxiter exit hands
    back whatever iterate it had, silently. Certify the stopping criterion
    on the returned solution instead (same rule scipy uses) and warn when it
    is not met — mirrors ``implicit_euler._check_solver_info`` on the numpy
    path, which raises on breakdown and warns on maxiter."""
    nb = float(jnp.linalg.norm(b))
    nr = float(jnp.linalg.norm(b - A(x)))
    if nr <= 10.0 * max(rtol * nb, atol):
        return
    logging.warning(
        "jax elliptic solve did not meet its stopping criterion: "
        f"|r|/|b| = {nr / max(nb, 1e-300):.2e} (rtol={rtol:g}, atol={atol:g}, "
        f"maxiter={maxiter}) — breakdown or maxiter; the iterate was used as is"
    )
