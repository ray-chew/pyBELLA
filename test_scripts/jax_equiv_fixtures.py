"""Realistic-state builders for the JAX operator-equivalence harness.

Numpy-only (no jax imports): builds small ModelStates from existing cases,
captures the exact linop-construction inputs the solver would pass
(reproducing the ``_prepare_{2d,3d}_system`` preamble up to, but excluding,
the scipy LinearOperator wrap), and provides deterministic test vectors.

Cases were picked to cover every static branch of the operators:
- travelling vortex 2D: periodic/periodic,
- internal long wave:   periodic/WALL (vertical wall zeroing),
- lamb wave:            ATMOSPHERIC_EXTENSION (y_atmosphere stencil wrap,
                        no slab zeroing in the divergence),
- smoke_agnesi 3D:      WALL vertical, with/without terrain metric
                        (contravariant fluxes, lap3D use_cross).
"""

import numpy as np

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import terrain
from pybella.flow_solver.numerics import coriolis, implicit_euler
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import (
    smoke_agnesi,
    test_internal_long_wave,
    test_lamb_wave,
    test_sphere_gw,
    test_sphere_swe_tc2,
)
from pybella.utils import axes, user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.operators.laplacian import preconditioner

import test_3d_elliptic_oracle as oracle3d


def _mem_from_case(case, **ud_overrides):
    ud = user_data.UserDataInit(**vars(case.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    for key, val in ud_overrides.items():
        setattr(ud, key, val)
    # production prepare.initialise calls this before grid_init (domain
    # extension for rayleigh sponge / forcing cases)
    if hasattr(ud, "rayleigh_bc"):
        ud.rayleigh_bc(ud)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = case.sol_init(sol, npf, elem, node, th, ud)
    # normally set by time_update from the eos switches; cases that don't
    # set it in sol_init get the unblended value
    if not hasattr(ud, "nonhydrostasy"):
        ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    if not hasattr(ud, "compressibility"):
        ud.compressibility = float(ud.is_compressible)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def make_vortex2d_mem():
    """2D travelling vortex, periodic/periodic, zero Coriolis/gravity."""
    ud = oracle3d.build_ud(iny=65, inz=1)
    mem = oracle3d.build_state(ud, "xy")
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def make_igw_mem():
    """Internal long wave: periodic x, WALL y."""
    return _mem_from_case(test_internal_long_wave)


def make_lamb_mem():
    """Lamb wave: ATMOSPHERIC_EXTENSION (y_atmosphere branch)."""
    return _mem_from_case(test_lamb_wave)


def make_agnesi3d_mem(with_terrain=True):
    """Quasi-2D 3D smoke_agnesi state, with or without its Agnesi hill."""
    overrides = {} if with_terrain else {"orography": None}
    mem, ud = _mem_from_case(smoke_agnesi, **overrides)
    assert (mem.elem.metric is not None) == with_terrain
    return mem, ud


def make_agnesi2d_mem():
    """Native-2D smoke_agnesi state with its Agnesi hill (inz=1)."""
    mem, ud = _mem_from_case(smoke_agnesi, inz=1)
    assert mem.elem.ndim == 2 and mem.elem.metric is not None
    return mem, ud


def make_sphere_swe_mem():
    """Thin-shell Williamson TC2: a non-vertical-line spherical metric with
    ``coriolis_field`` f(phi) e_r, general phi/r free-slip walls, e_up
    buoyancy and no gravity. Coarsened (pole-safe phi resolution) and with
    the initial projection skipped for fixture speed."""
    mem, ud = _mem_from_case(
        test_sphere_swe_tc2, inx=32 + 1, inz=48 + 1, initial_projection=False
    )
    assert mem.elem.metric is not None and not mem.elem.metric.vertical_line
    assert mem.elem.metric.e_up is not None
    assert getattr(ud, "coriolis_field", None) is not None
    return mem, ud


def make_sphere_gw_mem():
    """3D compressible spherical shell (DCMIP-31 gravity wave): the true
    r-dependent metric with radial gravity, e_up buoyancy (general H^-1),
    a general phi free-slip wall and the well-balanced gravity ghost fill on
    the radial axis. Coarsened for fixture speed."""
    mem, ud = _mem_from_case(test_sphere_gw, inx=32 + 1, iny=4 + 1, inz=8 + 1)
    assert mem.elem.metric is not None and not mem.elem.metric.vertical_line
    assert mem.elem.metric.e_up is not None
    return mem, ud


def capture_lap2d_args(mem, ud, dt=0.01):
    """Linop-construction args of ``_prepare_2d_system``: (npf, node,
    coriolis_params, diag_inv, ud)."""
    bdry_c.set_ghost_cells(mem, ud)
    implicit_euler.operator_coefficients_nodes(mem, ud, dt)

    coriolis_params = coriolis.multiply_inverse_terms(
        mem.npf, mem, ud, dt, attrs=("u", "v", "w"), get_coeffs=True
    )

    if mem.elem.metric is not None:
        h11_t, h22_t, h12_t, h21_t = coriolis_params
        h2x2 = ((h11_t.T, h12_t.T), (h21_t.T, h22_t.T))
        M = terrain.elliptic_tensor_2d(mem.elem.metric, h2x2)
        coriolis_params = (M[0][0].T, M[1][1].T, M[0][1].T, M[1][0].T)
        met = mem.elem.metric
        diag_inv = preconditioner.prepare_diag(
            mem.npf,
            mem.node,
            cii=(
                mem.npf.wplus[0] * met.J,
                mem.npf.wplus[1] * (1.0 + met.G1 * met.G1) * met.ooJ,
                None,
            ),
        )
    else:
        diag_inv = preconditioner.prepare_diag(mem.npf, mem.node)

    return (mem.npf, mem.node, coriolis_params, diag_inv, ud)


def capture_lap3d_args(mem, ud, dt=None):
    """Linop-construction args of ``_prepare_3d_system``: (elem, node, npf,
    ud, diag_inv, dt, cij)."""
    dt = float(ud.dtfixed) if dt is None else dt
    bdry_c.set_ghost_cells(mem, ud)
    implicit_euler.operator_coefficients_nodes(mem, ud, dt)

    hv = coriolis.compute_inverse_coefficients(mem, ud, dt)
    h_role = ((hv[0], hv[1], hv[2]), (hv[3], hv[4], hv[5]), (hv[6], hv[7], hv[8]))
    if mem.elem.metric is not None:
        h_role = terrain.elliptic_tensor(mem.elem.metric, h_role)
    rho_of = axes.role_of_axis(axes.vertical_axis(ud))
    cij = [
        [mem.npf.wplus[i] * h_role[rho_of[i]][rho_of[j]] for j in range(3)]
        for i in range(3)
    ]
    diag_inv = preconditioner.prepare_diag(
        mem.npf, mem.node, cii=(cij[0][0], cij[1][1], cij[2][2])
    )
    return (mem.elem, mem.node, mem.npf, ud, diag_inv, dt, cij)


def test_vectors(n, special_indices=(), seed=99):
    """Deterministic matvec probes: smooth, random, constant, plus delta
    spikes at the given flat indices (boundary-wrap bugs hide from smooth
    vectors)."""
    x = np.linspace(0.0, 2.0 * np.pi, n)
    vecs = [
        np.sin(3.0 * x) + 0.2 * np.cos(7.0 * x),
        np.random.default_rng(seed).standard_normal(n),
        np.ones(n),
    ]
    for i in special_indices:
        v = np.zeros(n)
        v[i] = 1.0
        vecs.append(v)
    return vecs


def lap2d_special_indices(node):
    """Corners, edge midpoints and centre of the row-major (iicy, iicx)
    solve grid — the nodes where x- and y-wraps compose."""
    nx, ny = node.iicx, node.iicy
    n = nx * ny
    return (
        0,
        nx - 1,
        n - nx,
        n - 1,
        nx // 2,
        n - nx + nx // 2,
        (ny // 2) * nx,
        (ny // 2) * nx + nx - 1,
        (ny // 2) * nx + nx // 2,
    )


def lap3d_special_indices(node):
    """Corners and centre of the C-ordered padded node.isc box."""
    sx, sy, sz = node.isc
    n = sx * sy * sz

    def flat(i, j, k):
        return (i * sy + j) * sz + k

    return (
        flat(1, 1, 1),
        flat(sx - 2, sy - 2, sz - 2),
        flat(1, sy - 2, 1),
        flat(sx - 2, 1, sz - 2),
        flat(sx // 2, sy // 2, sz // 2),
        0,
        n - 1,
    )
