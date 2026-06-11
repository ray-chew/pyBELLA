"""JAX backend equivalence: every jax_ops twin vs its numpy golden master.

Component 1 of the JAX migration (dev_notes/jax_migration.md): the pure
operator layer. Inputs are realistic coefficient fields captured from small
real cases (see jax_equiv_fixtures), covering periodic/WALL/atmosphere
boundaries and the terrain metric, plus synthetic fields for shape coverage.

Tolerances (x64), scaled by field magnitude: XLA's CPU codegen lowers
float64 division to reciprocal-multiply (verified: jax output equals
``a * (1/s)`` bitwise, and ``--xla_cpu_enable_fast_math=false`` does not
change it), so division-bearing ports deviate by ~1 ulp — ULP_TOL = 1e-14
relative. Ports whose reduction order can also differ get ACCUM = 1e-13.
Both are many orders tighter than the 1e-5..1e-8 the golden masters
themselves are held to, and float32 leakage (~1e-7) fails outright.

Skips cleanly when jax is not installed.
"""

import copy

import numpy as np
import pytest

jax = pytest.importorskip("jax")

from pybella.backends import jax_ops  # noqa: E402  (enables x64 at import)
from pybella.utils.operators import (  # noqa: E402
    convolution as conv_np,
    divergence as div_np,
    finite_difference as fd_np,
    gradient as grad_np,
)
from pybella.utils.operators.laplacian import (  # noqa: E402
    lap2D_manual as lap2D_np,
    lap3D as lap3D_np,
    preconditioner as precon_np,
)

import jax_equiv_fixtures as fx  # noqa: E402

ULP_TOL = 1e-14
ACCUM = 1e-13


def maxdiff(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))


def assert_close(got, want, tol, label=""):
    want = np.asarray(want)
    scale = max(1.0, float(np.max(np.abs(want))))
    diff = maxdiff(got, want)
    assert (
        diff <= tol * scale
    ), f"{label} max|diff| = {diff:.3e} > {tol:.0e} * {scale:.3e}"


def rand(shape, seed=0):
    return np.random.default_rng(seed).standard_normal(shape)


def test_x64_is_enabled():
    assert jax.config.jax_enable_x64
    assert jax.numpy.ones(1).dtype == np.float64


# ---------------------------------------------------------------- fixtures


@pytest.fixture(scope="session")
def vortex():
    return fx.make_vortex2d_mem()


@pytest.fixture(scope="session")
def igw():
    return fx.make_igw_mem()


@pytest.fixture(scope="session")
def lamb():
    return fx.make_lamb_mem()


@pytest.fixture(scope="session")
def agnesi_terrain():
    return fx.make_agnesi3d_mem(with_terrain=True)


@pytest.fixture(scope="session")
def agnesi_flat():
    return fx.make_agnesi3d_mem(with_terrain=False)


@pytest.fixture(scope="session")
def agnesi_2d():
    return fx.make_agnesi2d_mem()


# ------------------------------------------------------- finite_difference


@pytest.mark.parametrize("shape", [(12, 9), (8, 7, 6), (5, 31)])
def test_finite_difference(shape):
    f = rand(shape, seed=1)
    for axis in range(len(shape)):
        got = jax_ops.finite_difference.do_1d(jax.numpy.asarray(f), 0.37, axis=axis)
        want = fd_np.do_1d(f, 0.37, axis=axis)
        assert_close(got, want, ULP_TOL)


# ----------------------------------------------------------------- gradient


@pytest.mark.parametrize("shape", [(12, 9), (8, 7, 6)])
def test_gradient_plain(shape):
    f = rand(shape, seed=2)
    if len(shape) == 2:
        want = grad_np.compute_gradient_2d(f, 0.1, 0.2)
        got = jax_ops.gradient.compute_gradient_2d(jax.numpy.asarray(f), 0.1, 0.2)
    else:
        want = grad_np.compute_gradient_3d(f, 0.1, 0.2, 0.3)
        got = jax_ops.gradient.compute_gradient_3d(jax.numpy.asarray(f), 0.1, 0.2, 0.3)
    for g, w in zip(got, want):
        assert_close(g, w, ULP_TOL)


@pytest.mark.parametrize("shape", [(12, 9), (8, 7, 6)])
def test_gradient_at_nodes_synthetic(shape):
    f = rand(shape, seed=3)
    ndim = len(shape)
    want = grad_np.compute_at_nodes(f, ndim, (0.1, 0.2, 0.3))
    got = jax_ops.gradient.compute_at_nodes(jax.numpy.asarray(f), ndim, (0.1, 0.2, 0.3))
    for g, w in zip(got, want):
        assert_close(g, w, ULP_TOL)


def test_gradient_at_nodes_real(vortex):
    mem, ud = vortex
    p = mem.npf.p2_nodes
    want = grad_np.compute_at_nodes(p, mem.elem.ndim, mem.node.dxyz)
    got = jax_ops.gradient.compute_at_nodes(
        jax.numpy.asarray(p), mem.elem.ndim, mem.node.dxyz
    )
    for g, w in zip(got, want):
        assert_close(g, w, ULP_TOL)


# -------------------------------------------------------------- convolution


@pytest.mark.parametrize("ndim", [2, 3])
def test_flux_kernels_directional(ndim):
    shape = (14, 11) if ndim == 2 else (9, 8, 7)
    data = rand(shape, seed=4)
    for direction, kernel in conv_np.get_flux_kernels(ndim).items():
        want = conv_np.apply_directional_convolution(data, kernel, direction, ndim)
        got = jax_ops.convolution.apply_directional_convolution(
            jax.numpy.asarray(data), kernel, direction, ndim
        )
        assert_close(got, want, ACCUM, label=direction)


@pytest.mark.parametrize("ndim,width", [(2, 2), (2, 3), (3, 2), (3, 3)])
def test_averaging_kernel(ndim, width):
    shape = (14, 11) if ndim == 2 else (9, 8, 7)
    data = rand(shape, seed=5)
    kernel = conv_np.get_averaging_kernel(ndim, width=width)
    want = conv_np.apply_convolution_kernel(data, kernel)
    got = jax_ops.convolution.apply_convolution_kernel(jax.numpy.asarray(data), kernel)
    assert_close(got, want, ACCUM)


def test_convolution_real_coefficient_field(vortex):
    # the rhoY**cexp field of operator_coefficients_nodes
    mem, ud = vortex
    cexp = 2.0 - mem.th.gamm
    data = mem.sol.rhoY**cexp
    kernel = conv_np.get_averaging_kernel(2, width=2)
    want = conv_np.apply_convolution_kernel(data, kernel)
    got = jax_ops.convolution.apply_convolution_kernel(jax.numpy.asarray(data), kernel)
    assert_close(got, want, ACCUM)


# --------------------------------------------------------------- divergence


@pytest.mark.parametrize("shape", [(12, 9), (8, 7, 6)])
def test_divergence_pure_kernels(shape):
    if len(shape) == 2:
        u, v = rand(shape, 6), rand(shape, 7)
        want = div_np.compute_2d(u, v, 0.1, 0.2)
        got = jax_ops.divergence.compute_2d(
            jax.numpy.asarray(u), jax.numpy.asarray(v), 0.1, 0.2
        )
        assert_close(got, want, ULP_TOL)
    else:
        u, v, w = rand(shape, 6), rand(shape, 7), rand(shape, 8)
        want = div_np.compute_3d_sum(u, v, w, 0.1, 0.2, 0.3)
        got = jax_ops.divergence.compute_3d_sum(
            *(jax.numpy.asarray(f) for f in (u, v, w)), 0.1, 0.2, 0.3
        )
        assert_close(got, want, ULP_TOL)
        want_c = div_np.compute_3d_components(u, v, w, 0.1, 0.2, 0.3)
        got_c = jax_ops.divergence.compute_3d_components(
            *(jax.numpy.asarray(f) for f in (u, v, w)), 0.1, 0.2, 0.3
        )
        for g, w_ in zip(got_c, want_c):
            assert_close(g, w_, ULP_TOL)


def _divergence_case(mem, ud):
    """Run numpy compute_at_nodes on deep copies; jax twin on the originals."""
    elem = mem.elem
    ndim = elem.ndim

    sol_np = copy.deepcopy(mem.sol)
    rhs_np = np.zeros_like(mem.npf.rhs)
    rhs_np = div_np.compute_at_nodes(rhs_np, elem, sol_np, ud)

    if elem.metric is not None:
        m = elem.metric
        metric = (
            (m.J, m.G1, None, None, None)
            if ndim == 2
            else (m.J, m.G1, m.G2, m.vaxis, m.haxes)
        )
    else:
        metric = None

    momenta = (mem.sol.rhou, mem.sol.rhov, mem.sol.rhow if ndim == 3 else None)
    rhs_jx, momenta_jx = jax_ops.divergence.compute_at_nodes(
        tuple(None if f is None else jax.numpy.asarray(f) for f in momenta),
        jax.numpy.asarray(mem.sol.rho),
        jax.numpy.asarray(mem.sol.rhoY),
        ndim,
        (elem.dx, elem.dy, elem.dz),
        wall_dims=jax_ops.divergence.wall_zero_dims(ud, ndim),
        metric=metric,
    )

    assert_close(rhs_jx, rhs_np, ULP_TOL, label="rhs")
    for jx, np_field in zip(momenta_jx, (sol_np.rhou, sol_np.rhov, sol_np.rhow)):
        if jx is not None:
            assert_close(jx, np_field, ULP_TOL, label="momentum")


def test_divergence_at_nodes_vortex(vortex):
    _divergence_case(*vortex)


def test_divergence_at_nodes_igw_wall(igw):
    _divergence_case(*igw)


def test_divergence_at_nodes_lamb_atmosphere(lamb):
    _divergence_case(*lamb)


def test_divergence_at_nodes_terrain(agnesi_terrain):
    _divergence_case(*agnesi_terrain)


def test_divergence_at_nodes_terrain_2d(agnesi_2d):
    _divergence_case(*agnesi_2d)


def test_divergence_at_nodes_3d_flat(agnesi_flat):
    _divergence_case(*agnesi_flat)


# ----------------------------------------------------------- preconditioner


def _precon_case(mem, ud, cii=None):
    want = precon_np.prepare_diag(mem.npf, mem.node, cii=cii)
    got = jax_ops.laplacian.preconditioner.prepare_diag(mem.npf, mem.node, cii=cii)
    assert_close(got, want, ACCUM)


def test_preconditioner_2d(vortex):
    mem, ud = vortex
    fx.capture_lap2d_args(mem, ud)  # materialise wplus/wcenter
    _precon_case(mem, ud)


def test_preconditioner_3d_cii(agnesi_terrain):
    mem, ud = agnesi_terrain
    args = fx.capture_lap3d_args(mem, ud)
    cij = args[6]
    _precon_case(mem, ud, cii=(cij[0][0], cij[1][1], cij[2][2]))


# ------------------------------------------------------------------- lap3D


def _lap3d_case(mem, ud):
    args = fx.capture_lap3d_args(mem, ud)
    lap_np = lap3D_np.get_linop(*args)
    lap_jx = jax_ops.laplacian.lap3D.get_linop(*args)
    n = int(np.prod(mem.node.isc))
    nonzero = 0
    for k, v in enumerate(fx.test_vectors(n, fx.lap3d_special_indices(mem.node))):
        want = np.asarray(lap_np(v))
        got = np.asarray(lap_jx(v))
        assert_close(got, want, ACCUM, label=f"vector {k}")
        nonzero += np.max(np.abs(want)) > 0
    assert nonzero >= 5  # most probes must actually exercise the operator


def test_lap3d_terrain_cross(agnesi_terrain):
    mem, ud = agnesi_terrain
    _lap3d_case(mem, ud)


def test_lap3d_flat(agnesi_flat):
    mem, ud = agnesi_flat
    _lap3d_case(mem, ud)


# ------------------------------------------------------------------- lap2D


def _lap2d_case(mem, ud):
    args = fx.capture_lap2d_args(mem, ud)
    lap_np = lap2D_np.get_linop(*args)
    lap_jx = jax_ops.laplacian.lap2D.get_linop(*args)
    n = mem.node.iicx * mem.node.iicy
    for k, v in enumerate(fx.test_vectors(n, fx.lap2d_special_indices(mem.node))):
        want = np.asarray(lap_np(v))
        got = np.asarray(lap_jx(v))
        assert np.max(np.abs(want)) > 0, f"vector {k} does not exercise the operator"
        assert_close(got, want, ACCUM, label=f"vector {k}")


def test_lap2d_periodic_periodic(vortex):
    _lap2d_case(*vortex)


def test_lap2d_periodic_wall(igw):
    _lap2d_case(*igw)


def test_lap2d_atmosphere(lamb):
    _lap2d_case(*lamb)


def test_lap2d_terrain(agnesi_2d):
    _lap2d_case(*agnesi_2d)
