"""JAX elliptic solve: seam + solver equivalence vs scipy.

Two layers of validation:

1. **Residual certificate** (implementation-independent): on manufactured
   systems built from the captured production coefficients, the JAX
   bicgstab solution must satisfy scipy's stopping criterion
   ``||b - A x|| <= max(rtol ||b||, atol)`` measured with the *numpy*
   operator. Krylov iterates differ between implementations; the
   certificate is what "solved the same system" means.

2. **Seam equivalence** (end-to-end): ``implicit_euler.do_implicit_part``
   run with ``ud.backend = "jax"`` vs ``"numpy"`` on deep-copied states
   must agree on p2_nodes and all momenta to solver tolerance
   (magnitude-scaled 1e-5 — see SOLVER_TOL below for the measurement).

Coefficient assembly stays on the numpy path for both backends (see
backends/jax_ops/elliptic_solve.py), so the systems solved are
bit-identical; only the operator application and the Krylov loop differ.

Skips cleanly when jax is not installed.
"""

import copy

import numpy as np
import pytest

jax = pytest.importorskip("jax")

from pybella.backends import jax_ops  # noqa: E402  (enables x64 at import)
from pybella.flow_solver.numerics import implicit_euler  # noqa: E402
from pybella.utils.operators.laplacian import (  # noqa: E402
    lap2D_manual as lap2D_np,
    lap3D as lap3D_np,
)

import jax_equiv_fixtures as fx  # noqa: E402

RTOL = 1e-5  # scipy bicgstab default rtol, matched by the jax path
ATOL = 1e-8  # ud.tol
CERT_SAFETY = 5.0  # recursive-vs-true residual drift headroom
# Field agreement: when rtol*||b|| dominates the stopping criterion (e.g.
# the stratified terrain case), two converged solutions legitimately differ
# at ~rtol relative — measured: both backends at true ||r|| ~ 2.7e-5 against
# bound 2.9e-5, solutions 1.2e-5 apart relative. 1e-5 scaled also mirrors
# the CompareSol regression tolerance.
SOLVER_TOL = 1e-5


def _assert_fields_close(mem_a, mem_b, label):
    fields = [("npf", "p2_nodes"), ("sol", "rhou"), ("sol", "rhov"), ("sol", "rhow")]
    for holder, name in fields:
        a = getattr(getattr(mem_a, holder), name)
        b = getattr(getattr(mem_b, holder), name)
        scale = max(1.0, float(np.max(np.abs(a))))
        diff = float(np.max(np.abs(a - b)))
        assert (
            diff <= SOLVER_TOL * scale
        ), f"{label}/{name}: max|diff| = {diff:.3e} > {SOLVER_TOL:.0e} * {scale:.3e}"


def _run_both_backends(mem, ud, dt):
    mem_np = copy.deepcopy(mem)
    mem_jx = copy.deepcopy(mem)
    try:
        ud.backend = "numpy"
        implicit_euler.do_implicit_part(mem_np, ud, dt)
        ud.backend = "jax"
        implicit_euler.do_implicit_part(mem_jx, ud, dt)
    finally:
        ud.backend = "numpy"
    return mem_np, mem_jx


# ------------------------------------------------- residual certificate


def _certify(lap_np_apply, lap_jx, n, label):
    """JAX-solved x must satisfy scipy's criterion under the numpy operator."""
    x_exact = np.sin(np.linspace(0.0, 6.0 * np.pi, n)) + 0.1
    b = np.asarray(lap_np_apply(x_exact)).ravel()
    x = jax_ops.elliptic_solve.bicgstab(lap_jx, b, atol=ATOL, maxiter=6000)
    r = b - np.asarray(lap_np_apply(x)).ravel()
    rnorm, bnorm = np.linalg.norm(r), np.linalg.norm(b)
    bound = CERT_SAFETY * max(RTOL * bnorm, ATOL)
    assert rnorm <= bound, f"{label}: ||r|| = {rnorm:.3e} > {bound:.3e}"


def test_certificate_lap2d_periodic():
    mem, ud = fx.make_vortex2d_mem()
    args = fx.capture_lap2d_args(mem, ud)
    _certify(
        lap2D_np.get_linop(*args),
        jax_ops.laplacian.lap2D.get_linop(*args),
        mem.node.iicx * mem.node.iicy,
        "lap2D vortex",
    )


def test_certificate_lap3d_terrain():
    mem, ud = fx.make_agnesi3d_mem(True)
    args = fx.capture_lap3d_args(mem, ud)
    _certify(
        lap3D_np.get_linop(*args),
        jax_ops.laplacian.lap3D.get_linop(*args),
        int(np.prod(mem.node.isc)),
        "lap3D agnesi terrain",
    )


# ------------------------------------------------------ seam equivalence


def test_seam_vortex2d():
    mem, ud = fx.make_vortex2d_mem()
    mem_np, mem_jx = _run_both_backends(mem, ud, 0.01)
    _assert_fields_close(mem_np, mem_jx, "vortex2d")


def test_seam_igw_wall():
    mem, ud = fx.make_igw_mem()
    mem_np, mem_jx = _run_both_backends(mem, ud, float(ud.dtfixed))
    _assert_fields_close(mem_np, mem_jx, "igw")


def test_seam_lamb_atmosphere():
    mem, ud = fx.make_lamb_mem()
    mem_np, mem_jx = _run_both_backends(mem, ud, float(ud.dtfixed))
    _assert_fields_close(mem_np, mem_jx, "lamb")


def test_seam_agnesi3d_terrain():
    mem, ud = fx.make_agnesi3d_mem(True)
    mem_np, mem_jx = _run_both_backends(mem, ud, float(ud.dtfixed))
    _assert_fields_close(mem_np, mem_jx, "agnesi3d-terrain")


def test_seam_agnesi3d_flat():
    mem, ud = fx.make_agnesi3d_mem(False)
    mem_np, mem_jx = _run_both_backends(mem, ud, float(ud.dtfixed))
    _assert_fields_close(mem_np, mem_jx, "agnesi3d-flat")


def test_seam_agnesi2d_terrain():
    mem, ud = fx.make_agnesi2d_mem()
    mem_np, mem_jx = _run_both_backends(mem, ud, float(ud.dtfixed))
    _assert_fields_close(mem_np, mem_jx, "agnesi2d-terrain")


def test_seam_sphere_pole_global():
    """Pole-to-pole shell: the elliptic pole-ring collapse — the JAX lap3D
    one-sided pole rows wrapped in the Galerkin scatter/gather —
    reproduces the numpy collapse through a full do_implicit_part. The systems
    are bit-identical (assembly is numpy on both backends), so this agrees to
    solver tolerance and, being a single well-conditioned solve, to machine
    precision in practice."""
    mem, ud = fx.make_sphere_swe_global_mem()
    from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c

    bdry_c.set_ghost_cells(mem, ud)
    dt = float(ud.dtfixed) if getattr(ud, "dtfixed", 0.0) else 0.01
    mem_np, mem_jx = _run_both_backends(mem, ud, dt)
    assert mem_jx._pole_collapse is not None and mem_np._pole_collapse is not None
    _assert_fields_close(mem_np, mem_jx, "sphere-pole-global")


def test_numpy_backend_is_default():
    mem, ud = fx.make_vortex2d_mem()
    assert getattr(ud, "backend", None) == "numpy"
