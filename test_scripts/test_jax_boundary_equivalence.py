"""JAX boundary-fill equivalence (Phase A of the device-resident plan).

Layer 1 — micro-tests: the jnp slice-assign fills vs the raw np.pad modes
and custom callbacks (`_negative_symmetric`, `periodic_plus_one`) on random
arrays; the custom callbacks have subtle reversed-slice semantics and
bit-faithfulness is the risk.

Layer 2 — per-BC equivalence: the public entry points (`set_ghost_cells`,
`set_ghost_nodes`, `rayleigh_damping`, `scale_wall_node_values`) run with
``ud.backend`` toggled on deep-copied states, parametrized over the fixture
cases (periodic / WALL+gravity / RAYLEIGH+ATMOSPHERIC_EXTENSION / terrain
3D + native-2D), step=None vs flipped sweep orientation, and the sol=
override. Pure explicit arithmetic — magnitude-scaled 1e-13.

Skips cleanly when jax is not installed.
"""

import copy

import numpy as np
import pytest

jax = pytest.importorskip("jax")

from pybella.backends import jax_ops  # noqa: E402  (enables x64 at import)
from pybella.flow_solver.utils.boundary import (  # noqa: E402
    cell_boundary as bdry_c,
    common as bdry_common,
    node_boundary as bdry_n,
    rayleigh_boundary as bdry_r,
)
from pybella.utils import axes  # noqa: E402

import jax_equiv_fixtures as fx  # noqa: E402

TOL = 1e-13


def assert_close(got, want, tol=TOL, label=""):
    want = np.asarray(want)
    scale = max(1.0, float(np.max(np.abs(want))))
    diff = float(np.max(np.abs(np.asarray(got) - want)))
    assert (
        diff <= tol * scale
    ), f"{label}: max|diff| = {diff:.3e} > {tol:.0e} * {scale:.3e}"


def rand(shape, seed):
    return np.random.default_rng(seed).standard_normal(shape)


# ------------------------------------------------------------- micro-tests


@pytest.mark.parametrize(
    "shape,dim", [((9, 7), 0), ((9, 7), 1), ((6, 5, 8), 2), ((8, 7, 5), 2)]
)
def test_micro_wrap_symmetric_negsym(shape, dim):
    from pybella.backends.jax_ops.boundary import (
        _negative_symmetric_fill,
        _symmetric_fill,
        _wrap_fill,
    )

    ig = 2
    f = rand(shape, seed=dim)
    pads = [(0, 0)] * len(shape)
    pads[dim] = (ig, ig)
    inner = [slice(None)] * len(shape)
    inner[dim] = slice(ig, -ig)
    inner = tuple(inner)

    want = np.pad(f[inner], pads, "wrap")
    assert_close(_wrap_fill(jax.numpy.asarray(f), dim, ig), want, label="wrap")

    want = np.pad(f[inner], pads, "symmetric")
    assert_close(
        _symmetric_fill(jax.numpy.asarray(f), dim, ig), want, label="symmetric"
    )

    want = np.pad(f[inner], pads, bdry_c._negative_symmetric)
    assert_close(
        _negative_symmetric_fill(jax.numpy.asarray(f), dim, ig), want, label="negsym"
    )


@pytest.mark.parametrize(
    "shape,dim", [((11, 8), 0), ((11, 8), 1), ((7, 9, 8), 1), ((7, 6, 9), 1)]
)
def test_micro_node_pads(shape, dim):
    from pybella.backends.jax_ops.boundary import (
        _periodic_plus_one_fill,
        _reflect_fill,
    )

    ig = 2
    f = rand(shape, seed=10 + dim)
    pads = [(0, 0)] * len(shape)
    pads[dim] = (ig, ig)
    inner = [slice(None)] * len(shape)
    inner[dim] = slice(ig, -ig)
    inner = tuple(inner)

    want = np.pad(f[inner], pads, bdry_n.periodic_plus_one)
    assert_close(
        _periodic_plus_one_fill(jax.numpy.asarray(f), dim, ig),
        want,
        label="periodic_plus_one",
    )

    want = np.pad(f[inner], pads, "reflect")
    assert_close(_reflect_fill(jax.numpy.asarray(f), dim, ig), want, label="reflect")


# ----------------------------------------------------- set_ghost_cells twins


def _scramble_ghosts(mem, seed=4):
    """Overwrite ghost regions with noise so fills must actually work."""
    rng = np.random.default_rng(seed)
    igs = mem.elem.igs
    for name in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX"):
        f = getattr(mem.sol, name)
        for dim in range(f.ndim):
            lo = [slice(None)] * f.ndim
            hi = [slice(None)] * f.ndim
            lo[dim] = slice(0, igs[dim])
            hi[dim] = slice(-igs[dim], None)
            f[tuple(lo)] = 0.5 + 0.1 * rng.standard_normal(f[tuple(lo)].shape)
            f[tuple(hi)] = 0.5 + 0.1 * rng.standard_normal(f[tuple(hi)].shape)
    # rho-derived fields must stay physical (positive rho/rhoY)
    np.abs(mem.sol.rho, out=mem.sol.rho)
    np.abs(mem.sol.rhoY, out=mem.sol.rhoY)


def _cells_case(mem, ud, step=None, flip=0, sol_override=False):
    mem_np = copy.deepcopy(mem)
    mem_jx = copy.deepcopy(mem)
    for m in (mem_np, mem_jx):
        _scramble_ghosts(m)
        for _ in range(flip):
            m.sol.flip_forward()
            if m.elem.metric is not None:
                m.elem.metric.flip_forward()

    def run(m, backend):
        try:
            ud.backend = backend
            sol = copy.deepcopy(m.sol) if sol_override else None
            bdry_c.set_ghost_cells(m, ud, step=step, sol=sol)
            return sol if sol_override else m.sol
        finally:
            ud.backend = "numpy"

    sol_np = run(mem_np, "numpy")
    sol_jx = run(mem_jx, "jax")
    for name in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX"):
        assert_close(
            getattr(sol_jx, name), getattr(sol_np, name), label=f"cells/{name}"
        )


def test_cells_vortex_periodic():
    _cells_case(*fx.make_vortex2d_mem())


def test_cells_igw_gravity_wall():
    _cells_case(*fx.make_igw_mem())


def test_cells_lamb_atmosphere_rayleigh():
    _cells_case(*fx.make_lamb_mem())


def test_cells_agnesi3d_terrain():
    _cells_case(*fx.make_agnesi3d_mem(True))


def test_cells_agnesi3d_flat():
    _cells_case(*fx.make_agnesi3d_mem(False))


def test_cells_agnesi2d_terrain():
    _cells_case(*fx.make_agnesi2d_mem())


def test_cells_sol_override():
    _cells_case(*fx.make_igw_mem(), sol_override=True)


@pytest.mark.parametrize("split", [0, 1])
def test_cells_sweep_flipped_igw(split):
    # mimic the advection sweep: split+1 forward flips, fill last axis
    mem, ud = fx.make_igw_mem()
    _cells_case(mem, ud, step=split, flip=split + 1)


@pytest.mark.parametrize("split", [0, 1, 2])
def test_cells_sweep_flipped_terrain3d(split):
    mem, ud = fx.make_agnesi3d_mem(True)
    _cells_case(mem, ud, step=split, flip=split + 1)


@pytest.mark.parametrize("split", [0, 1])
def test_cells_sweep_flipped_lamb(split):
    mem, ud = fx.make_lamb_mem()
    _cells_case(mem, ud, step=split, flip=split + 1)


def test_cells_incompressible_branch():
    mem, ud = fx.make_igw_mem()
    ud.is_compressible = 0
    _cells_case(mem, ud)


# ----------------------------------------------- general (spherical) metric


def test_cells_sphere_swe_walls():
    # thin-shell TC2: general free-slip mirror on the phi + degenerate-r walls
    _cells_case(*fx.make_sphere_swe_mem())


def test_cells_sphere_gw_gravity_and_wall():
    # 3D shell: well-balanced e_up gravity fill on r + general phi wall mirror
    _cells_case(*fx.make_sphere_gw_mem())


def test_cells_sphere_swe_sweep_phi():
    # phi wall reached mid-sweep at the extremal split (canonical orientation)
    mem, ud = fx.make_sphere_swe_mem()
    _cells_case(mem, ud, step=2, flip=3)


def test_cells_sphere_gw_sweep_phi():
    mem, ud = fx.make_sphere_gw_mem()
    _cells_case(mem, ud, step=2, flip=3)


def test_cells_sphere_gw_sweep_gravity():
    # radial (gravity) sweep: general e_up fill in the sweep orientation
    mem, ud = fx.make_sphere_gw_mem()
    _cells_case(mem, ud, step=1, flip=2)


# ------------------------------------------------------- lat-lon pole fold


def test_cells_sphere_pole_global():
    # full pole-to-pole TC2: the phi POLE fold (pure index remap of all 6
    # fields) + the periodic lambda + degenerate-r general wall, canonical
    _cells_case(*fx.make_sphere_swe_global_mem())


def test_cells_sphere_pole_sweep_phi():
    # the phi POLE fold reached mid-sweep (phi swept last: split=2, 3 flips)
    mem, ud = fx.make_sphere_swe_global_mem()
    _cells_case(mem, ud, step=2, flip=3)


# ------------------------------------------------------ set_ghost_nodes twin


def _nodes_case(mem, ud):
    p_np = rand(mem.npf.p2_nodes.shape, seed=7)
    p_jx = p_np.copy()
    try:
        ud.backend = "numpy"
        bdry_n.set_ghost_nodes(p_np, mem.node, ud)
        ud.backend = "jax"
        bdry_n.set_ghost_nodes(p_jx, mem.node, ud)
    finally:
        ud.backend = "numpy"
    assert_close(p_jx, p_np, label="nodes")


def test_nodes_vortex():
    _nodes_case(*fx.make_vortex2d_mem())


def test_nodes_igw():
    _nodes_case(*fx.make_igw_mem())


def test_nodes_agnesi3d_quasi2d():
    _nodes_case(*fx.make_agnesi3d_mem(True))


def test_nodes_sphere_gw():
    _nodes_case(*fx.make_sphere_gw_mem())


def test_nodes_sphere_pole_global():
    # pole node fold: gather (reflect about the pole node) + pole-row ring-mean
    _nodes_case(*fx.make_sphere_swe_global_mem())


# ------------------------------------------------- rayleigh + scale-wall twins


def test_rayleigh_damping_sponge():
    mem, ud = fx.make_lamb_mem()
    mem_np, mem_jx = copy.deepcopy(mem), copy.deepcopy(mem)
    try:
        ud.backend = "numpy"
        bdry_r.rayleigh_damping(mem_np.sol, mem_np.npf, ud)
        ud.backend = "jax"
        bdry_r.rayleigh_damping(mem_jx.sol, mem_jx.npf, ud)
    finally:
        ud.backend = "numpy"
    for name in ("rhou", "rhov", "rhow", "rhoY"):
        assert_close(
            getattr(mem_jx.sol, name), getattr(mem_np.sol, name), label=f"ray/{name}"
        )


def test_rayleigh_damping_with_forcing_arrays():
    from pybella.tests import test_unstable_lamb

    mem, ud = fx._mem_from_case(test_unstable_lamb)
    assert hasattr(ud, "rayleigh_forcing") and ud.rayleigh_forcing
    # synthetic forcing arrays with the right shapes; func-mode semantics
    shp = mem.sol.rho.shape
    forcing = [
        rand(shp, 1),
        rand(shp, 2),
        rand(shp, 3),
        rand(mem.npf.p2_nodes.shape, 4),
        0.3,
    ]
    mem_np, mem_jx = copy.deepcopy(mem), copy.deepcopy(mem)
    try:
        ud.backend = "numpy"
        bdry_r.rayleigh_damping(mem_np.sol, mem_np.npf, ud, forcing=list(forcing))
        ud.backend = "jax"
        bdry_r.rayleigh_damping(mem_jx.sol, mem_jx.npf, ud, forcing=list(forcing))
    finally:
        ud.backend = "numpy"
    for name in ("rhou", "rhov", "rhow", "rhoY"):
        assert_close(
            getattr(mem_jx.sol, name), getattr(mem_np.sol, name), label=f"rayf/{name}"
        )
    assert_close(mem_jx.npf.p2_nodes, mem_np.npf.p2_nodes, label="rayf/p2_nodes")


def test_scale_wall_node_values():
    mem, ud = fx.make_igw_mem()
    rhs_np = rand(mem.npf.rhs.shape, seed=9)
    rhs_jx = rhs_np.copy()
    try:
        ud.backend = "numpy"
        bdry_common.scale_wall_node_values(rhs_np, mem.node, ud)
        ud.backend = "jax"
        bdry_common.scale_wall_node_values(rhs_jx, mem.node, ud)
    finally:
        ud.backend = "numpy"
    assert_close(rhs_jx, rhs_np, label="scale_wall")


# --------------------------------------------------------------- spy test


def test_jax_fill_path_is_live():
    from unittest import mock

    mem, ud = fx.make_igw_mem()
    calls = {"n": 0}
    orig = jax_ops.boundary.set_ghost_cells

    def spy(*a, **k):
        calls["n"] += 1
        return orig(*a, **k)

    try:
        ud.backend = "jax"
        with mock.patch.object(jax_ops.boundary, "set_ghost_cells", side_effect=spy):
            bdry_c.set_ghost_cells(mem, ud)
    finally:
        ud.backend = "numpy"
    assert calls["n"] == 1
