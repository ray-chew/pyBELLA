"""JAX advection / coriolis / diffusion equivalence (component 3).

The backend seam for advection sits inside
``compute_advection._compute_flux_and_recovery``: the numpy driver (sweeps,
flips, ghost fills, flux-difference updates) is shared, only the per-sweep
recovery+HLL kernel and the advective-flux convolution are swapped. These
tests run the *public drivers* with ``ud.backend`` toggled on deep-copied
states and compare every solution field — pure explicit arithmetic, so
tolerances are at the XLA-ulp/accumulation floor (magnitude-scaled), far
below solver tolerances.

Skips cleanly when jax is not installed.
"""

import copy

import numpy as np
import pytest

jax = pytest.importorskip("jax")

from pybella.backends import jax_ops  # noqa: E402  (enables x64 at import)
from pybella.flow_solver.numerics import coriolis as coriolis_np  # noqa: E402
from pybella.flow_solver.numerics import diffusion as diffusion_np  # noqa: E402
from pybella.flow_solver.numerics.explicit_advection import (  # noqa: E402
    advective_flux,
    compute_advection,
)

import jax_equiv_fixtures as fx  # noqa: E402

TOL = 1e-12  # several chained ulp-level kernels per sweep


def assert_close(got, want, tol, label=""):
    want = np.asarray(want)
    scale = max(1.0, float(np.max(np.abs(want))))
    diff = float(np.max(np.abs(np.asarray(got) - want)))
    assert (
        diff <= tol * scale
    ), f"{label}: max|diff| = {diff:.3e} > {tol:.0e} * {scale:.3e}"


def _assert_sol_close(sol_jx, sol_np, label):
    for name in ("rho", "rhou", "rhov", "rhow", "rhoX", "rhoY"):
        assert_close(
            getattr(sol_jx, name), getattr(sol_np, name), TOL, f"{label}/{name}"
        )


def _run_both(mem, ud, fn):
    """Run fn(mem, ud) once per backend on deep copies; return both mems."""
    mem_np = copy.deepcopy(mem)
    mem_jx = copy.deepcopy(mem)
    try:
        ud.backend = "numpy"
        fn(mem_np, ud)
        ud.backend = "jax"
        fn(mem_jx, ud)
    finally:
        ud.backend = "numpy"
    return mem_np, mem_jx


# -------------------------------------------------- advective flux (rhoY)


def _flux_case(mem, ud):
    def fn(m, u):
        advective_flux.recompute(m, u)

    mem_np, mem_jx = _run_both(mem, ud, fn)
    fl_np = mem_np.cache.get_flux_containers(mem_np.elem)
    fl_jx = mem_jx.cache.get_flux_containers(mem_jx.elem)
    for dim in range(mem.elem.ndim):
        assert_close(fl_jx[dim].rhoY, fl_np[dim].rhoY, TOL, f"flux[{dim}].rhoY")


def test_advective_flux_vortex():
    _flux_case(*fx.make_vortex2d_mem())


def test_advective_flux_terrain():
    _flux_case(*fx.make_agnesi3d_mem(True))


# ------------------------------------------------------- advection drivers


def _strang_case(mem, ud, dt):
    def fn(m, u):
        advective_flux.recompute(m, u)
        compute_advection.strange_splitting(m, u, dt, odd=0, label="equiv")
        compute_advection.strange_splitting(m, u, dt, odd=1, label="equiv")

    mem_np, mem_jx = _run_both(mem, ud, fn)
    _assert_sol_close(mem_jx.sol, mem_np.sol, "strang")


def _rk_case(mem, ud, dt):
    def fn(m, u):
        advective_flux.recompute(m, u)
        compute_advection.first_order_runge_kutta(m, u, 0.5 * dt)

    mem_np, mem_jx = _run_both(mem, ud, fn)
    _assert_sol_close(mem_jx.sol, mem_np.sol, "rk")
    assert_close(mem_jx.sol.pwchi, mem_np.sol.pwchi, TOL, "rk/pwchi")


def test_strang_vortex2d():
    mem, ud = fx.make_vortex2d_mem()
    _strang_case(mem, ud, 0.01)


def test_strang_igw_wall():
    mem, ud = fx.make_igw_mem()
    _strang_case(mem, ud, float(ud.dtfixed))


def test_strang_agnesi3d_terrain():
    mem, ud = fx.make_agnesi3d_mem(True)
    _strang_case(mem, ud, float(ud.dtfixed))


def test_strang_agnesi2d_terrain():
    mem, ud = fx.make_agnesi2d_mem()
    _strang_case(mem, ud, float(ud.dtfixed))


def test_rk_vortex2d():
    mem, ud = fx.make_vortex2d_mem()
    _rk_case(mem, ud, 0.01)


def test_rk_agnesi3d_terrain():
    mem, ud = fx.make_agnesi3d_mem(True)
    _rk_case(mem, ud, float(ud.dtfixed))


# ---------------------------------------------------------------- coriolis


def test_coriolis_multiply_inverse():
    mem, ud = fx.make_igw_mem()
    ud.coriolis_strength = np.array([0.2, 0.3, 0.5])
    dt = float(ud.dtfixed)

    def fn(m, u):
        coriolis_np.multiply_inverse_terms(m.sol, m, u, dt)

    mem_np, mem_jx = _run_both(mem, ud, fn)
    for name in ("rhou", "rhov", "rhow"):
        assert_close(getattr(mem_jx.sol, name), getattr(mem_np.sol, name), 1e-13, name)


def test_coriolis_get_coeffs_and_inverse_coefficients():
    mem, ud = fx.make_igw_mem()
    ud.coriolis_strength = np.array([0.0, 0.0, 0.4])
    dt = float(ud.dtfixed)

    mem_np = copy.deepcopy(mem)
    mem_jx = copy.deepcopy(mem)
    try:
        ud.backend = "numpy"
        c_np = coriolis_np.multiply_inverse_terms(
            mem_np.npf, mem_np, ud, dt, attrs=("u", "v", "w"), get_coeffs=True
        )
        v_np = coriolis_np.compute_inverse_coefficients(mem_np, ud, dt)
        ud.backend = "jax"
        c_jx = coriolis_np.multiply_inverse_terms(
            mem_jx.npf, mem_jx, ud, dt, attrs=("u", "v", "w"), get_coeffs=True
        )
        v_jx = coriolis_np.compute_inverse_coefficients(mem_jx, ud, dt)
    finally:
        ud.backend = "numpy"

    for k, (a, b) in enumerate(zip(c_jx, c_np)):
        assert_close(a, b, 1e-13, f"coeff2d[{k}]")
    for k, (a, b) in enumerate(zip(v_jx, v_np)):
        assert_close(a, b, 1e-13, f"hinv[{k}]")


def _sphere_coriolis_case(mem, ud, label):
    """Compare apply + inverse coefficients on a non-vertical-line metric
    (general H^-1 kernel): the buoyancy rank-one term along e_up and, for
    TC2, the spatially varying ``coriolis_field`` role components."""
    dt = float(getattr(ud, "dtfixed", 5.0)) or 5.0

    mem_np = copy.deepcopy(mem)
    mem_jx = copy.deepcopy(mem)
    try:
        ud.backend = "numpy"
        coriolis_np.multiply_inverse_terms(mem_np.sol, mem_np, ud, dt)
        v_np = coriolis_np.compute_inverse_coefficients(mem_np, ud, dt)
        v_np = [np.asarray(x).copy() for x in v_np]
        ud.backend = "jax"
        coriolis_np.multiply_inverse_terms(mem_jx.sol, mem_jx, ud, dt)
        v_jx = coriolis_np.compute_inverse_coefficients(mem_jx, ud, dt)
        v_jx = [np.asarray(x).copy() for x in v_jx]
    finally:
        ud.backend = "numpy"

    for name in ("rhou", "rhov", "rhow"):
        assert_close(
            getattr(mem_jx.sol, name),
            getattr(mem_np.sol, name),
            1e-13,
            f"{label}/apply/{name}",
        )
    for k, (a, b) in enumerate(zip(v_jx, v_np)):
        assert_close(a, b, 1e-13, f"{label}/hinv[{k}]")


def test_coriolis_general_sphere_swe():
    # thin-shell TC2: coriolis_field f(phi) e_r + e_up buoyancy, general H^-1
    _sphere_coriolis_case(*fx.make_sphere_swe_mem(), "sphere_swe")


def test_coriolis_general_sphere_gw():
    # 3D shell gravity wave: e_up buoyancy (nu != 0), zero rotation
    _sphere_coriolis_case(*fx.make_sphere_gw_mem(), "sphere_gw")


# --------------------------------------------------------------- diffusion


def _diffusion_case(mem, ud):
    ud.diffusion = True
    ud.diffusion_coeff = 1.0e-4
    dt = float(getattr(ud, "dtfixed", 0.01)) or 0.01

    def fn(m, u):
        diffusion_np.apply(m, u, dt)

    mem_np, mem_jx = _run_both(mem, ud, fn)
    _assert_sol_close(mem_jx.sol, mem_np.sol, "diffusion")


def test_diffusion_2d():
    _diffusion_case(*fx.make_igw_mem())


def test_diffusion_3d():
    _diffusion_case(*fx.make_agnesi3d_mem(False))
