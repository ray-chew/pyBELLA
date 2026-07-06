"""Stage-A gates for the spherical lat-lon metric (sphere pt 0).

The sphere enters as metric data from ``SphericalShellMap`` (a Tier-3
``CurvilinearMap``: vertical coordinate lines are radial, not parallel
Cartesian lines). These gates pin, on the spherical channel:

1. closed forms: J = r~^2 cos(phi), N_r = J e_r, e_up = e_r, |t_r| = 1,
   height = r - a, G1/G2 unset, coordinates = the documented embedding;
2. duality N_a . t_b = J delta_ab and det N = J^2;
3. the discrete metric identity sum_a D_a N_a: -> 0 at 2nd order for the
   true map; for the frozen-radius (thin-shell/SWE) map the defect is
   PURELY RADIAL and equals -2 a cos(phi) e_r — tangential fluxes are
   never contaminated (the D1 shell-degeneracy contract);
4. gradient-map completeness: apply_gradient_map inverts the exact
   chain rule to roundoff (sum_a (N_a)_k (t_a)_j = J delta_kj);
5. flux-divergence oracle: discrete sum_a D_a(N_a . f) -> J div f at
   2nd order for an analytic Cartesian vector field;
6. elliptic fold on the sphere: M = (1/J) N N^T symmetric positive
   definite, and the assembled elliptic operator == div o correction
   (the map-agnostic composition identity of
   ``test_terrain_elliptic_oracle.py``, rerun on the spherical metric);
7. activation: ud.curvilinear_map dispatches in grid_init; combining
   with ud.orography raises; the frozen-radius metric is exactly
   r-uniform (the quasi-2D broadcast prerequisite); flips round-trip
   ``e_up``.
"""

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import spherical, terrain
from pybella.utils import options as opts

# ---------------------------------------------------------------- helpers

_A = 2.0  # nondimensional planet radius for the metric-level gates


class _GridStubUD:
    """Minimal ud for grid_init on the spherical channel."""

    def __init__(self, n=(32, 8, 16), a=_A, frozen=False, phi_max=1.2, cmap=True):
        self.inx, self.iny, self.inz = n[0] + 1, n[1] + 1, n[2] + 1
        self.xmin, self.xmax = -np.pi, np.pi  # lambda
        self.ymin, self.ymax = a * 0.95, a * 1.05  # r
        self.zmin, self.zmax = -phi_max, phi_max  # phi
        self.bdry_type = np.array(
            [opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.WALL]
        )
        self.gravity_direction = 1
        if cmap:
            self.curvilinear_map = spherical.SphericalShellMap(a, frozen_radius=frozen)


def _coord(grid_obj, axis):
    from pybella.utils import axes as _axes

    shape = [1] * grid_obj.ndim
    shape[axis] = -1
    return _axes.coords_along(grid_obj, axis).reshape(shape)


def _grid_coords(grid_obj):
    """(lam, r, phi) broadcast views on one grid."""
    return tuple(_coord(grid_obj, a) for a in range(3))


def _e_r(lam, phi):
    """Radial unit vector of the documented embedding (pole along -x_2)."""
    cp = np.cos(phi)
    return (cp * np.cos(lam), cp * np.sin(lam), -np.sin(phi) + 0.0 * lam)


def _metrics(n=(32, 8, 16), frozen=False, phi_max=1.2):
    # ghost latitudes must stay inside (-pi/2, pi/2); coarse convergence
    # grids pass a smaller phi_max so the J > 0 guard is not (correctly)
    # tripped by ghost rows beyond the pole
    ud = _GridStubUD(n=n, frozen=frozen, phi_max=phi_max)
    elem, node = dis_grid.grid_init(ud)
    return ud, elem, node


# --------------------------------------------------- closed forms + duality


@pytest.mark.parametrize("frozen", [False, True], ids=["true", "frozen"])
@pytest.mark.parametrize("loc", [0, 1], ids=["cells", "nodes"])
def test_spherical_closed_forms(frozen, loc):
    ud, elem, node = _metrics(frozen=frozen)
    grid_obj = (elem, node)[loc]
    m = grid_obj.metric
    lam, r, phi = _grid_coords(grid_obj)
    r_eff = _A if frozen else r

    assert m.vertical_line is False
    assert m.G1 is None and m.G2 is None

    J_expect = np.broadcast_to(r_eff**2 * np.cos(phi), m.J.shape)
    np.testing.assert_allclose(m.J, J_expect, rtol=1e-14)
    assert np.all(m.J > 0.0)

    # N_r = J e_r and e_up = e_r
    er = _e_r(lam, phi)
    for k in range(3):
        np.testing.assert_allclose(
            m.N[1][k],
            np.broadcast_to(m.J * er[k], m.J.shape),
            rtol=1e-13,
            atol=1e-14 * _A**2,
        )
        np.testing.assert_allclose(
            m.e_up[k], np.broadcast_to(er[k], m.J.shape), rtol=1e-13, atol=1e-15
        )

    # |t_r| = 1 (unit radial tangent) and height = r - a
    np.testing.assert_allclose(m.h_v, np.ones_like(m.J), rtol=1e-14)
    np.testing.assert_allclose(
        m.height, np.broadcast_to(r - _A, m.J.shape), rtol=0, atol=1e-14 * _A
    )

    # coordinates: the documented embedding, all three materialized
    cp = np.cos(phi)
    x_expect = (r * cp * np.cos(lam), r * cp * np.sin(lam), -r * np.sin(phi))
    for k in range(3):
        assert m.x[k] is not None
        np.testing.assert_allclose(
            m.x[k], np.broadcast_to(x_expect[k], m.J.shape), rtol=1e-14, atol=1e-15
        )


@pytest.mark.parametrize("frozen", [False, True], ids=["true", "frozen"])
def test_spherical_duality_and_detN(frozen):
    ud, elem, _ = _metrics(frozen=frozen)
    m = elem.metric
    t = [
        [np.broadcast_to(c, m.J.shape) for c in ta]
        for ta in ud.curvilinear_map.tangents(_grid_coords(elem))
    ]
    scale = np.max(np.abs(m.J))
    for a in range(3):
        for b in range(3):
            expect = m.J if a == b else 0.0
            np.testing.assert_allclose(
                sum(m.N[a][k] * t[b][k] for k in range(3)),
                expect,
                rtol=1e-13,
                atol=1e-14 * scale,
            )
    # det N = J^2 (nonsingularity)
    N = m.N
    detN = sum(
        N[0][k]
        * (
            N[1][(k + 1) % 3] * N[2][(k + 2) % 3]
            - N[1][(k + 2) % 3] * N[2][(k + 1) % 3]
        )
        for k in range(3)
    )
    np.testing.assert_allclose(detN, m.J * m.J, rtol=1e-12)


# ------------------------------------------- discrete metric identity


def _div_normals(elem):
    """Interior central differences sum_a D_a N_a, per Cartesian k."""
    m = elem.metric
    d = (elem.dx, elem.dy, elem.dz)
    inner = (slice(1, -1),) * 3
    out = []
    for k in range(3):
        acc = np.zeros_like(m.J[inner])
        for a in range(3):
            arr = m.N[a][k]
            sl_p = [slice(1, -1)] * 3
            sl_m = [slice(1, -1)] * 3
            sl_p[a] = slice(2, None)
            sl_m[a] = slice(0, -2)
            acc += (arr[tuple(sl_p)] - arr[tuple(sl_m)]) / (2.0 * d[a])
        out.append(acc)
    return out


def test_metric_identity_true_map_second_order():
    errs = []
    for n in ((16, 4, 8), (32, 8, 16)):
        _, elem, _ = _metrics(n=n, frozen=False, phi_max=1.0)
        div = _div_normals(elem)
        errs.append(max(np.max(np.abs(c)) for c in div) / _A**2)
    assert errs[1] < errs[0] / 3.0  # ~4 for 2nd order
    assert errs[1] < 2e-2


def test_metric_identity_frozen_defect_is_radial():
    """Frozen shell: sum_a D_a N_a = -2 a cos(phi) e_r + O(h^2); the
    tangential projection converges to ZERO (no tangential-flux bias)."""
    tang_errs, rad_errs = [], []
    for n in ((16, 4, 8), (32, 8, 16)):
        _, elem, _ = _metrics(n=n, frozen=True, phi_max=1.0)
        div = _div_normals(elem)
        inner = (slice(1, -1),) * 3
        lam, r, phi = _grid_coords(elem)
        lam_i = np.broadcast_to(lam, elem.metric.J.shape)[inner]
        phi_i = np.broadcast_to(phi, elem.metric.J.shape)[inner]
        er = _e_r(lam_i, phi_i)
        radial = sum(div[k] * er[k] for k in range(3))
        tang = [div[k] - radial * er[k] for k in range(3)]
        tang_errs.append(max(np.max(np.abs(c)) for c in tang) / _A)
        rad_errs.append(np.max(np.abs(radial - (-2.0 * _A * np.cos(phi_i)))) / _A)
    assert tang_errs[1] < tang_errs[0] / 3.0
    assert rad_errs[1] < rad_errs[0] / 3.0
    assert tang_errs[1] < 5e-3 and rad_errs[1] < 5e-2


# --------------------------------------------- gradient map + divergence


@pytest.mark.parametrize("frozen", [False, True], ids=["true", "frozen"])
def test_gradient_map_inverts_chain_rule(frozen):
    """With EXACT computational gradients dp/dxi_a = t_a . grad_x p, the
    gradient map must return grad_x p to roundoff — the completeness
    identity sum_a (N_a)_k (t_a)_j = J delta_kj, algebraic (no h)."""
    ud, elem, _ = _metrics(frozen=frozen)
    m = elem.metric
    x = m.x
    # p = x0^2 + sin(x1) + x2 x0 -> grad = (2 x0 + x2, cos(x1), x0)
    grad_cart = (2.0 * x[0] + x[2], np.cos(x[1]), x[0])
    t = [
        [np.broadcast_to(c, m.J.shape) for c in ta]
        for ta in ud.curvilinear_map.tangents(_grid_coords(elem))
    ]
    dp = [sum(t[a][k] * grad_cart[k] for k in range(3)) for a in range(3)]
    out = terrain.apply_gradient_map(m, list(dp))
    scale = max(np.max(np.abs(g)) for g in grad_cart)
    for k in range(3):
        np.testing.assert_allclose(out[k], grad_cart[k], rtol=1e-12, atol=1e-13 * scale)


def test_flux_divergence_oracle_second_order():
    """Discrete sum_a D_a (N_a . f) -> J div f for analytic f(x)."""
    errs = []
    for n in ((16, 4, 8), (32, 8, 16)):
        _, elem, _ = _metrics(n=n, frozen=False, phi_max=1.0)
        m = elem.metric
        x = m.x
        f = (x[0] ** 2, np.sin(x[1]), x[2] * x[0])
        divf = 3.0 * x[0] + np.cos(x[1])

        d = (elem.dx, elem.dy, elem.dz)
        inner = (slice(1, -1),) * 3
        acc = np.zeros_like(m.J[inner])
        for a in range(3):
            F = sum(m.N[a][k] * f[k] for k in range(3))
            sl_p = [slice(1, -1)] * 3
            sl_m = [slice(1, -1)] * 3
            sl_p[a] = slice(2, None)
            sl_m[a] = slice(0, -2)
            acc += (F[tuple(sl_p)] - F[tuple(sl_m)]) / (2.0 * d[a])
        expect = (m.J * divf)[inner]
        errs.append(np.max(np.abs(acc - expect)) / np.max(np.abs(expect)))
    assert errs[1] < errs[0] / 3.0
    assert errs[1] < 5e-2


# ------------------------------------------------------- elliptic fold


def test_elliptic_fold_spd_on_sphere():
    _, elem, _ = _metrics(frozen=False)
    m = elem.metric
    one = np.ones_like(m.J)
    zero = np.zeros_like(m.J)
    ident = ((one, zero, zero), (zero, one, zero), (zero, zero, one))
    M = terrain.elliptic_tensor(m, ident)

    scale = np.max(np.abs(M[0][0]))
    for r in range(3):
        for s in range(r + 1, 3):
            np.testing.assert_allclose(M[r][s], M[s][r], rtol=1e-13, atol=1e-14 * scale)

    d1 = M[0][0]
    d2 = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    d3 = (
        M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
        - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
        + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
    )
    assert np.all(d1 > 0.0) and np.all(d2 > 0.0) and np.all(d3 > 0.0)


def test_elliptic_composition_identity_on_sphere():
    """The assembled elliptic operator == div o correction on the spherical
    channel (zero gravity: Stage A2 owns the gravity-boundary path)."""
    from pybella.flow_solver.numerics import implicit_euler
    from pybella.flow_solver.physics import hydrostatics, thermodynamics
    from pybella.flow_solver.utils import cache, fields
    from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
    from pybella.flow_solver.utils.boundary import node_boundary as bdry_n
    from pybella.tests import smoke_agnesi
    from pybella.utils import axes, user_data
    from pybella.utils.data_structures import ModelState
    from pybella.utils.operators import divergence
    from pybella.utils.operators.laplacian import preconditioner
    from pybella.flow_solver.numerics import coriolis

    base = smoke_agnesi.UserData()
    base.grav = 0.0
    base.gravity_strength = np.zeros(3)
    base.stratification = lambda y: 1.0 + 0.0 * y
    base.orography = None
    base.xmin, base.xmax = -np.pi, np.pi
    base.ymin, base.ymax = _A * 0.95, _A * 1.05
    base.zmin, base.zmax = -1.2, 1.2
    base.bdry_type[2] = opts.BdryType.WALL
    base.inx, base.iny, base.inz = 32 + 1, 8 + 1, 16 + 1
    base.curvilinear_map = spherical.SphericalShellMap(_A, frozen_radius=False)
    ud = user_data.UserDataInit(**vars(base))
    ud.coriolis_strength = np.array(ud.coriolis_strength)

    elem, node = dis_grid.grid_init(ud)
    assert elem.metric is not None and elem.metric.vertical_line is False

    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    hydrostatics.integrated_state(npf, elem, node, th, ud)

    # smooth deterministic state (any consistent fields do; g = 0)
    m = elem.metric
    sol.rho[...] = 1.0 + 0.1 * np.sin(m.x[0]) * np.cos(m.x[2])
    sol.rhoY[...] = 1.0 + 0.05 * np.cos(m.x[1])
    sol.rhou[...] = 0.1 * np.sin(m.x[2])
    sol.rhov[...] = 0.05 * np.cos(m.x[0])
    sol.rhow[...] = 0.02 * np.sin(m.x[1])
    sol.rhoX[...] = 0.0
    npf.p2_nodes[...] = 0.0
    ud.nonhydrostasy = 1.0
    ud.compressibility = 1.0

    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)  # zero gravity: no G1/G2 path

    dt = float(ud.dtfixed)
    implicit_euler.operator_coefficients_nodes(mem, ud, dt)
    lap, _ = implicit_euler._prepare_linear_system(mem, ud, dt)

    x = node.x[1:-1].reshape(-1, 1, 1)
    y = node.y[1:-1].reshape(1, -1, 1)
    z = node.z[1:-1].reshape(1, 1, -1)
    Lx = node.x[-1] - node.x[0]
    p_box = (
        np.sin(2 * np.pi * x / Lx) * np.cos(np.pi * y)
        + 0.1 * np.cos(2 * np.pi * x / Lx + 0.3) * y
        + 0.05 * np.sin(z)
    )
    lhs = np.asarray(lap @ p_box.ravel()).reshape(node.isc)

    hv = coriolis.compute_inverse_coefficients(mem, ud, dt)
    h_role = ((hv[0], hv[1], hv[2]), (hv[3], hv[4], hv[5]), (hv[6], hv[7], hv[8]))
    h_role = terrain.elliptic_tensor(mem.elem.metric, h_role)
    rho_of = axes.role_of_axis(axes.vertical_axis(ud))
    cij = [
        [mem.npf.wplus[i] * h_role[rho_of[i]][rho_of[j]] for j in range(3)]
        for i in range(3)
    ]
    diag_inv = preconditioner.prepare_diag(
        mem.npf, mem.node, cii=(cij[0][0], cij[1][1], cij[2][2])
    )

    p_full = np.zeros(node.sc)
    p_full[node.i1] = p_box
    bdry_n.set_ghost_nodes(p_full, node, ud)

    mem.sol.rhou[...] = 0.0
    mem.sol.rhov[...] = 0.0
    mem.sol.rhow[...] = 0.0
    implicit_euler._correction_nodes(mem, ud, dt, p_full, 0)
    rhs = np.zeros(node.isc)
    divergence.compute_at_nodes(rhs, mem.elem, mem.sol, ud)

    comp = diag_inv * (-(1.0 / dt) * rhs + mem.npf.wcenter * p_box)

    win = tuple(slice(3, -3) if s > 8 else slice(1, -1) for s in lhs.shape)
    scale = np.max(np.abs(lhs[win]))
    assert np.max(np.abs((lhs - comp)[win])) / scale <= 1e-12


# ----------------------------------------------- activation + mechanics


def test_activation_dispatch_and_guard():
    # no map, no orography -> metric is None (bit-identity contract)
    ud_flat = _GridStubUD(cmap=False)
    elem, node = dis_grid.grid_init(ud_flat)
    assert elem.metric is None and node.metric is None

    # map + orography -> hard error until SphericalTerrainMap composes them
    ud_bad = _GridStubUD()
    ud_bad.orography = lambda xi1, xi2: 0.0 * xi1
    with pytest.raises(ValueError, match="curvilinear_map"):
        dis_grid.grid_init(ud_bad)


def test_frozen_radius_metric_is_r_uniform():
    """The thin-shell degeneracy: every frozen-shell metric array is
    exactly uniform along the radial array axis (broadcast-safe)."""
    _, elem, node = _metrics(frozen=True)
    for grid_obj in (elem, node):
        m = grid_obj.metric
        arrays = [m.J, m.h_v] + [c for Na in m.N for c in Na] + list(m.e_up)
        for arr in arrays:
            assert np.array_equal(arr, np.broadcast_to(arr[:, :1, :], arr.shape))


def test_flip_cycle_preserves_e_up():
    _, elem, _ = _metrics(frozen=False)
    m = elem.metric
    e0 = [c.copy() for c in m.e_up]
    m.flip_forward()
    m.flip_backward()
    for k in range(3):
        assert np.array_equal(m.e_up[k], e0[k])
    for _ in range(3):
        m.flip_forward()
    for k in range(3):
        assert np.array_equal(m.e_up[k], e0[k])
