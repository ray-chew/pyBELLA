"""Reduction contract for the general curvilinear metric.

Executable definition of "reduces to the vertical-line terrain code"
for the general metric (J, N1, N2, N3): for the vertical-line map
z = z(xi_h1, eta, xi_h2) (x, y identity) the tangents in role order
(h1, v, h2) are

    t1 = (1, G1, 0),   t2 = (0, J, 0),   t3 = (0, G2, 1)

with J = dz/deta and G_h = dz/dxi_h, and the area normals N_i = t_j x t_k
(cyclic) must come out as

    N1 = (J, 0, 0),   N2 = (-G1, 1, -G2),   N3 = (0, 0, J)

so the general flux components F_i = N_i . f reproduce
``_metric_contravariant_fluxes_jit`` BIT-EXACTLY:

    F1 = J f_h1,   F2 = f_v - G1 f_h1 - G2 f_h2,   F3 = J f_h2.

The bit-exact contraction order is part of the contract: vertical
component first, then h1, then h2 (see ``general_flux``) — any general
flux assembly must keep this order so the reduction stays exact.

Also asserted: the duality N_i . t_j = J delta_ij (equivalent to
N_i = J grad xi_i, the bridge to the gradient/elliptic picture),
det N = J^2 (nonsingularity), and a horizontally stretched fixture
x = (x(xi1), z(xi1, eta, xi3), y(xi3)) pinning the identities the
stretched-grid advection and wall-reflection tests depend on — in
particular the effective-slope identity -(N2)_h / (N2)_v = z_xi_h / x'_h
that keeps the contravariant wall reflection
(cell_boundary._slope_terms) correct on stretched grids.
"""

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import terrain
from pybella.tests import smoke_agnesi
from pybella.utils import user_data
from pybella.utils.operators import divergence

# ---------------------------------------------------------------- helpers


def cross(a, b):
    """Elementwise cross product of role-component triples of arrays."""
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def normals(t1, t2, t3):
    """Area normals N_i = t_j x t_k for cyclic (i, j, k)."""
    return cross(t2, t3), cross(t3, t1), cross(t1, t2)


def general_flux(N_i, f):
    """The contraction contract: F_i = N_i . f, vertical first.

    Role-component order is (h1, v, h2); summing v + h1 + h2 left-to-right
    makes the reduction bit-exact against the vertical-line kernel
    ``divergence._metric_contravariant_fluxes_jit`` (1*f_v == f_v,
    x + (-y) == x - y, 0 + x == x in IEEE arithmetic).
    """
    return N_i[1] * f[1] + N_i[0] * f[0] + N_i[2] * f[2]


def vertical_line_tangents(m):
    """Role-component tangents of the current vertical-only map."""
    one = np.ones_like(m.J)
    zero = np.zeros_like(m.J)
    G2 = m.G2 if m.G2 is not None else zero
    t1 = (one, m.G1, zero)
    t2 = (zero, m.J, zero)
    t3 = (zero, G2, one)
    return t1, t2, t3


def _agnesi_hill(ud, h0_m=300.0, a_m=5000.0):
    h0 = h0_m / ud.h_ref
    a = a_m / ud.h_ref

    def h(xi1, xi2):
        # genuinely 2D orography so G1 and G2 are both nontrivial
        return h0 * a**2 / (xi1**2 + a**2) * (1.0 + 0.3 * np.sin(3.0 * xi2))

    return h


def _make_metrics(sleve=False, native_2d=False):
    """(elem.metric, node.metric) on the smoke_agnesi grid with a hill."""
    ud = user_data.UserDataInit(**vars(smoke_agnesi.UserData()))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    if native_2d:
        ud.inz = 1
    ud.orography = _agnesi_hill(ud)
    if sleve:
        hill = ud.orography
        ud.orography_smooth = lambda xi1, xi2: 0.5 * hill(xi1, xi2)
        ud.vertical_transform = terrain.SLEVETransform(
            s1=6000.0 / ud.h_ref, s2=1500.0 / ud.h_ref
        )
    elem, node = dis_grid.grid_init(ud)
    return elem.metric, node.metric


def _random_fields(rng, shape):
    rho = 1.0 + 0.5 * rng.random(shape)
    rhoY = 0.8 + 0.4 * rng.random(shape)
    moms = tuple(rng.random(shape) - 0.5 for _ in range(3))
    return rho, rhoY, moms


# ------------------------------------------------- vertical-line maps


@pytest.mark.parametrize("sleve", [False, True], ids=["galchen", "sleve"])
@pytest.mark.parametrize("loc", [0, 1], ids=["cells", "nodes"])
def test_vertical_line_normals_and_duality(sleve, loc):
    m = _make_metrics(sleve=sleve)[loc]
    t1, t2, t3 = vertical_line_tangents(m)
    N1, N2, N3 = normals(t1, t2, t3)

    # J reduction: J_full = t1 . (t2 x t3) == J, exactly
    assert np.array_equal(dot(t1, N1), m.J)

    # the three normals reduce to the documented forms, exactly
    zero = np.zeros_like(m.J)
    one = np.ones_like(m.J)
    assert all(np.array_equal(a, b) for a, b in zip(N1, (m.J, zero, zero)))
    assert all(np.array_equal(a, b) for a, b in zip(N2, (-m.G1, one, -m.G2)))
    assert all(np.array_equal(a, b) for a, b in zip(N3, (zero, zero, m.J)))

    # duality N_i . t_j = J delta_ij, exactly (G1 - G1 == 0 etc.)
    T = (t1, t2, t3)
    N = (N1, N2, N3)
    for i in range(3):
        for j in range(3):
            expect = m.J if i == j else zero
            assert np.array_equal(dot(N[i], T[j]), expect)

    # nonsingularity: det of the normal matrix is J^2, exactly
    detN = dot(N1, cross(N2, N3))
    assert np.array_equal(detN, m.J * m.J)


@pytest.mark.parametrize("sleve", [False, True], ids=["galchen", "sleve"])
def test_flux_reduction_3d_bit_exact(sleve):
    """General F_i = N_i . (theta m) == today's kernel, bit for bit."""
    m = _make_metrics(sleve=sleve)[0]
    N1, N2, N3 = normals(*vertical_line_tangents(m))

    rng = np.random.default_rng(20260612)
    rho, rhoY, moms = _random_fields(rng, m.J.shape)
    mom_h1, mom_v, mom_h2 = moms

    f_h1_leg, f_v_leg, f_h2_leg = divergence._metric_contravariant_fluxes_jit(
        rho, rhoY, mom_h1, mom_v, mom_h2, m.J, m.G1, m.G2
    )

    theta = rhoY / rho
    f = (mom_h1 * theta, mom_v * theta, mom_h2 * theta)
    assert np.array_equal(general_flux(N1, f), f_h1_leg)
    assert np.array_equal(general_flux(N2, f), f_v_leg)
    assert np.array_equal(general_flux(N3, f), f_h2_leg)


def test_flux_reduction_2d_bit_exact():
    """2D restriction: N1 = (J, 0), N2 = (-G1, 1) against the 2D kernel."""
    m = _make_metrics(native_2d=True)[0]
    assert m.G2 is None and m.J.ndim == 2

    # 2D area normals: N_i = J grad xi_i = perp of the other tangent
    t1 = (np.ones_like(m.J), m.G1)
    t2 = (np.zeros_like(m.J), m.J)
    N1 = (t2[1], -t2[0])
    N2 = (-t1[1], t1[0])

    # reduction + duality (J delta_ij), exactly
    assert np.array_equal(N1[0], m.J) and np.array_equal(N2[0], -m.G1)
    det2 = lambda a, b: a[0] * b[1] - a[1] * b[0]
    assert np.array_equal(det2(t1, t2), m.J)  # J = det of tangents
    assert np.array_equal(N1[0] * t1[0] + N1[1] * t1[1], m.J)
    assert np.array_equal(N2[0] * t2[0] + N2[1] * t2[1], m.J)
    assert np.array_equal(N1[0] * t2[0] + N1[1] * t2[1], np.zeros_like(m.J))
    assert np.array_equal(N2[0] * t1[0] + N2[1] * t1[1], np.zeros_like(m.J))

    rng = np.random.default_rng(20260612)
    rho, rhoY, (mom_h1, mom_v, _) = _random_fields(rng, m.J.shape)
    f_h1_leg, f_v_leg = divergence._metric_contravariant_fluxes_2d_jit(
        rho, rhoY, mom_h1, mom_v, m.J, m.G1
    )

    theta = rhoY / rho
    f = (mom_h1 * theta, mom_v * theta)
    assert np.array_equal(N2[1] * f[1] + N2[0] * f[0], f_v_leg)
    assert np.array_equal(N1[1] * f[1] + N1[0] * f[0], f_h1_leg)


# ------------------------------------- horizontally stretched fixture


def _tier2_map():
    """Analytic stretched map x = (x(xi1), z(xi1, eta, xi3), y(xi3)).

    Gal-Chen z over a 2D hill, with smooth monotone horizontal stretches;
    returns role-component tangents plus the analytic pieces.
    """
    xi1 = np.linspace(-1.5, 1.5, 24)
    eta = np.linspace(0.0, 1.0, 20)
    xi3 = np.linspace(-0.8, 0.8, 16)
    X1, E, X3 = np.meshgrid(xi1, eta, xi3, indexing="ij")

    xp = 1.0 + 0.35 * np.cos(1.3 * X1)  # x'(xi1) > 0
    yp = 1.0 + 0.25 * np.sin(0.9 * X3)  # y'(xi3) > 0

    h = 0.1 / (1.0 + X1**2) * (1.0 + 0.3 * np.cos(2.0 * X3))
    h_x1 = -0.2 * X1 / (1.0 + X1**2) ** 2 * (1.0 + 0.3 * np.cos(2.0 * X3))
    h_x3 = 0.1 / (1.0 + X1**2) * (-0.6 * np.sin(2.0 * X3))

    b = 1.0 - E  # Gal-Chen decay, eta0 = 0, etat = 1
    z_eta = 1.0 - h  # dz/deta (eta-independent for Gal-Chen)
    z_x1 = h_x1 * b
    z_x3 = h_x3 * b

    zero = np.zeros_like(E)
    t1 = (xp, z_x1, zero)
    t2 = (zero, z_eta, zero)
    t3 = (zero, z_x3, yp)
    return t1, t2, t3, xp, yp, z_eta, z_x1, z_x3


def test_tier2_tangent_fixture():
    t1, t2, t3, xp, yp, z_eta, z_x1, z_x3 = _tier2_map()
    N1, N2, N3 = normals(t1, t2, t3)
    J = dot(t1, N1)

    # closed forms of this map's normals and Jacobian
    np.testing.assert_allclose(J, xp * z_eta * yp, rtol=1e-14)
    assert np.all(J > 0.0)
    np.testing.assert_allclose(N1[0], z_eta * yp, rtol=1e-14)
    assert np.array_equal(N1[1], np.zeros_like(J))
    np.testing.assert_allclose(N2[0], -(yp * z_x1), rtol=1e-14)
    np.testing.assert_allclose(N2[1], xp * yp, rtol=1e-14)
    np.testing.assert_allclose(N2[2], -(xp * z_x3), rtol=1e-14)
    np.testing.assert_allclose(N3[2], xp * z_eta, rtol=1e-14)

    # duality and nonsingularity on the stretched map (real FP cancellation)
    T = (t1, t2, t3)
    N = (N1, N2, N3)
    scale = np.max(np.abs(J))
    for i in range(3):
        for j in range(3):
            expect = J if i == j else np.zeros_like(J)
            np.testing.assert_allclose(
                dot(N[i], T[j]), expect, rtol=1e-13, atol=1e-14 * scale
            )
    np.testing.assert_allclose(dot(N1, cross(N2, N3)), J * J, rtol=1e-13)

    # effective-slope identity: -(N2)_h / (N2)_v = z_xi_h / x'_h — the
    # quantity cell_boundary's contravariant wall reflection relies on
    np.testing.assert_allclose(-N2[0] / N2[1], z_x1 / xp, rtol=1e-13)
    np.testing.assert_allclose(-N2[2] / N2[1], z_x3 / yp, rtol=1e-13)

    # reflection invariance: N2 . m / (N2)_v is the effective-slope
    # contravariant momentum (positive per-column factor divides out)
    rng = np.random.default_rng(42)
    m = tuple(rng.random(J.shape) - 0.5 for _ in range(3))
    contra_eff = m[1] - (z_x1 / xp) * m[0] - (z_x3 / yp) * m[2]
    np.testing.assert_allclose(dot(N2, m) / N2[1], contra_eff, rtol=1e-12, atol=1e-13)


# ------------------------------------------- MetricFields N/x machinery
#
# Optional: these exercise the generalized metric (terrain.CurvilinearMap).
# They skip cleanly on a tree that only has the algebraic contract above.

_general_metric = pytest.mark.skipif(
    not hasattr(terrain, "CurvilinearMap"),
    reason="general curvilinear-map metric machinery not present",
)


@_general_metric
@pytest.mark.parametrize("sleve", [False, True], ids=["galchen", "sleve"])
@pytest.mark.parametrize("loc", [0, 1], ids=["cells", "nodes"])
def test_builder_normals_match_cross_products(sleve, loc):
    """MetricFields.N (synthesized) == cross products of the tangents,
    bit for bit (smoke_agnesi is v = 1, so role order == Cartesian order)."""
    m = _make_metrics(sleve=sleve)[loc]
    N_expect = normals(*vertical_line_tangents(m))
    a_h1, a_h2 = m.haxes
    for axis, Ni in zip((a_h1, m.vaxis, a_h2), N_expect):
        for k in range(3):
            assert np.array_equal(m.N[axis][k], Ni[k])
    # x: only the vertical coordinate is materialized (identity elsewhere)
    assert m.x[m.cart_v] is m.z
    assert m.x[a_h1] is None and m.x[a_h2] is None
    assert (m.cart_v, m.cart_haxes) == (m.vaxis, m.haxes)


@_general_metric
def test_flip_rotates_normals_consistently():
    """Sweep flips roll leaf axes and rotate WHICH normal sits on which
    array axis; Cartesian components never permute, so the vertical
    normal's slope components track the (rolled) G1/G2 slope arrays."""
    m = _make_metrics()[0]
    ndim = m.J.ndim
    N0 = [[c.copy() for c in Na] for Na in m.N]
    vax0, hax0 = m.vaxis, m.haxes

    m.flip_forward()
    one = np.ones_like(m.J)
    # the vertical normal now lives at the shifted vaxis; components are
    # still Cartesian: cart_h1 slot carries -G1 (rolled with the leaves)
    assert np.array_equal(m.N[m.vaxis][m.cart_haxes[0]], -m.G1)
    assert np.array_equal(m.N[m.vaxis][m.cart_v], one)
    assert np.array_equal(m.N[m.vaxis][m.cart_haxes[1]], -m.G2)
    assert np.array_equal(m.N[m.haxes[0]][m.cart_haxes[0]], m.J)
    assert np.array_equal(m.x[m.cart_v], m.z)

    # full cycle of flips is the identity (exact)
    for _ in range(ndim - 1):
        m.flip_forward()
    assert (m.vaxis, m.haxes) == (vax0, hax0)
    for a in range(ndim):
        for k in range(ndim):
            assert np.array_equal(m.N[a][k], N0[a][k])

    # and backward inverts forward
    m.flip_forward()
    m.flip_backward()
    for a in range(ndim):
        for k in range(ndim):
            assert np.array_equal(m.N[a][k], N0[a][k])


# --------------------------------------- general path (CurvilinearMap)


class _StubUD:
    """Minimal ud for grid_init: all-WALL so coordinate wraps are identity
    (required for the bit-exact comparison against build_metric_fields)."""

    def __init__(self, orography=None, orography_grad=None):
        self.inx, self.iny, self.inz = 17, 9, 7
        self.xmin, self.xmax = -1.0, 1.0
        self.ymin, self.ymax = 0.0, 2.0
        self.zmin, self.zmax = -0.5, 0.5
        from pybella.utils import options as opts

        self.bdry_type = np.array([opts.BdryType.WALL] * 3)
        self.gravity_direction = 1
        if orography is not None:
            self.orography = orography
        if orography_grad is not None:
            self.orography_grad = orography_grad


def _stub_hill(amplitude=0.05, width=0.3):
    def h(xi1, xi2):
        return amplitude / (1.0 + (xi1 / width) ** 2) * (1.0 + 0.2 * xi2)

    def dh1(xi1, xi2):
        return (
            -2.0 * amplitude * xi1 / width**2 / (1.0 + (xi1 / width) ** 2) ** 2
        ) * (1.0 + 0.2 * xi2)

    def dh2(xi1, xi2):
        return amplitude / (1.0 + (xi1 / width) ** 2) * 0.2 + 0.0 * xi1

    return h, (dh1, dh2)


def _coord(grid_obj, axis):
    from pybella.utils import axes as _axes

    shape = [1] * grid_obj.ndim
    shape[axis] = -1
    return _axes.coords_along(grid_obj, axis).reshape(shape)


if hasattr(terrain, "CurvilinearMap"):

    class _GalChenLineMap(terrain.CurvilinearMap):
        """The vertical-line Gal-Chen map written as a CurvilinearMap, using
        the same closed forms as ``build_metric_fields`` (bit-exact)."""

        def __init__(self, hill, grads, eta0, etat):
            self.h, (self.dh1, self.dh2) = hill, grads
            self.eta0, self.etat = eta0, etat

        def _decay(self, eta):
            return (self.etat - eta) / (self.etat - self.eta0)

        def coordinates(self, xi):
            xi1, eta, xi2 = xi
            z = eta + self.h(xi1, xi2) * self._decay(eta)
            return [None, z, None]

        def tangents(self, xi):
            xi1, eta, xi2 = xi
            b = self._decay(eta)
            J = (1.0 - self.h(xi1, xi2) / (self.etat - self.eta0)) + 0.0 * eta
            G1 = self.dh1(xi1, xi2) * b
            G2 = self.dh2(xi1, xi2) * b
            zero = 0.0 * J
            one = 1.0 + zero
            return [[one, G1, zero], [zero, J, zero], [zero, G2, one]]

    class _Tier2StretchMap(_GalChenLineMap):
        """x-stretched Gal-Chen: x = s(xi1), z = eta + h b(eta), y = xi3."""

        def __init__(self, hill, grads, eta0, etat, ax=0.3, kx=1.1):
            super().__init__(hill, grads, eta0, etat)
            self.ax, self.kx = ax, kx

        def coordinates(self, xi):
            xi1, eta, xi2 = xi
            _, z, _ = super().coordinates(xi)
            x = xi1 + self.ax * np.sin(self.kx * xi1)
            return [x + 0.0 * eta + 0.0 * xi2, z, None]

        def tangents(self, xi):
            xi1, eta, xi2 = xi
            t = super().tangents(xi)
            xp = 1.0 + self.ax * self.kx * np.cos(self.kx * xi1)
            t[0][0] = xp + 0.0 * t[0][1]
            return t


@_general_metric
def test_general_path_reduces_to_legacy_builder():
    """build_metric_fields_from_map == build_metric_fields, bit for bit,
    for the vertical-line Gal-Chen map."""
    hill, grads = _stub_hill()
    ud = _StubUD(orography=hill, orography_grad=grads)
    elem, node = dis_grid.grid_init(ud)
    cmap = _GalChenLineMap(hill, grads, eta0=ud.ymin, etat=ud.ymax)
    for grid_obj in (elem, node):
        legacy = grid_obj.metric
        general = terrain.build_metric_fields_from_map(grid_obj, ud, cmap)
        assert np.array_equal(general.J, legacy.J)
        assert np.array_equal(general.ooJ, legacy.ooJ)
        assert np.array_equal(general.G1, legacy.G1)
        assert np.array_equal(general.G2, legacy.G2)
        assert np.array_equal(general.z, legacy.z)
        for a in range(3):
            for k in range(3):
                assert np.array_equal(general.N[a][k], legacy.N[a][k])
        assert (general.vaxis, general.haxes) == (legacy.vaxis, legacy.haxes)


@_general_metric
def test_general_path_tier2_effective_slopes():
    """On an x-stretched grid the general builder must produce J = x'.J_z,
    duality-consistent normals, and effective slopes G1_eff = z_xi1 / x'."""
    hill, grads = _stub_hill()
    ud = _StubUD(orography=hill, orography_grad=grads)
    elem, _ = dis_grid.grid_init(ud)
    cmap = _Tier2StretchMap(hill, grads, eta0=ud.ymin, etat=ud.ymax)
    m = terrain.build_metric_fields_from_map(elem, ud, cmap)

    t = [
        [np.broadcast_to(c, m.J.shape) for c in ta]
        for ta in cmap.tangents([_coord(elem, a) for a in range(3)])
    ]
    # duality N_a . t_b = J delta_ab on the genuinely stretched map
    scale = np.max(np.abs(m.J))
    for a in range(3):
        for b_ in range(3):
            expect = m.J if a == b_ else 0.0
            np.testing.assert_allclose(
                sum(m.N[a][k] * t[b_][k] for k in range(3)),
                expect,
                rtol=1e-13,
                atol=1e-14 * scale,
            )
    # effective slopes: -(N_v)_h / (N_v)_v == z_xi_h / x'_h
    np.testing.assert_allclose(m.G1, t[0][1] / t[0][0], rtol=1e-13, atol=1e-15)
    # J = x' * dz/deta * 1
    np.testing.assert_allclose(m.J, t[0][0] * t[1][1], rtol=1e-13)


@_general_metric
def test_elliptic_fold_spd_symmetry():
    """M = (1/J) N H^-1 N^T with symmetric H^-1 = I must be symmetric and
    positive definite (det N = J^2 != 0) on a genuinely stretched map —
    the solvability condition of the projection."""
    hill, grads = _stub_hill()
    ud = _StubUD(orography=hill, orography_grad=grads)
    elem, _ = dis_grid.grid_init(ud)
    cmap = _Tier2StretchMap(hill, grads, eta0=ud.ymin, etat=ud.ymax)
    m = terrain.build_metric_fields_from_map(elem, ud, cmap)

    one = np.ones_like(m.J)
    zero = np.zeros_like(m.J)
    ident = ((one, zero, zero), (zero, one, zero), (zero, zero, one))
    M = terrain.elliptic_tensor(m, ident)

    scale = np.max(np.abs(M[0][0]))
    for r in range(3):
        for s in range(r + 1, 3):
            np.testing.assert_allclose(M[r][s], M[s][r], rtol=1e-13, atol=1e-14 * scale)

    # Sylvester minors pointwise: M is PD wherever J > 0
    d1 = M[0][0]
    d2 = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    d3 = (
        M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
        - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
        + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
    )
    assert np.all(d1 > 0.0) and np.all(d2 > 0.0) and np.all(d3 > 0.0)


@_general_metric
def test_advective_flux_upwind_sign_consistency():
    """Flux and upwind decision are monotone in the same signed quantity:
    F_i = (rhoY/rho)(N_i . m) and the Courant velocity F_i/(rhoY J) carry
    sign(N_i . v) — no sign inconsistency possible.
    Random fields, both signs, on a genuinely stretched map."""
    from pybella.flow_solver.numerics.explicit_advection import advective_flux

    hill, grads = _stub_hill()
    ud = _StubUD(orography=hill, orography_grad=grads)
    elem, _ = dis_grid.grid_init(ud)
    cmap = _Tier2StretchMap(hill, grads, eta0=ud.ymin, etat=ud.ymax)
    m = terrain.build_metric_fields_from_map(elem, ud, cmap)

    rng = np.random.default_rng(7)
    shape = m.J.shape

    class _Sol:
        pass

    sol = _Sol()
    sol.rho = 1.0 + 0.5 * rng.random(shape)
    sol.rhoY = 0.8 + 0.4 * rng.random(shape)
    sol.rhou, sol.rhov, sol.rhow = (rng.random(shape) - 0.5 for _ in range(3))

    for i in range(3):
        contra = advective_flux._normal_momentum(sol, m, i)
        flux = sol.rhoY * contra / sol.rho
        courant = flux / (sol.rhoY * m.J)
        # all three share the sign of N_i . m (rho, rhoY, J > 0)
        assert np.array_equal(np.sign(flux), np.sign(contra))
        assert np.array_equal(np.sign(courant), np.sign(contra))
        # and the Courant velocity is the contravariant speed (N_i . v)/J
        np.testing.assert_allclose(
            courant, contra / (sol.rho * m.J), rtol=1e-13, atol=1e-16
        )


@_general_metric
def test_general_path_rejects_nonpositive_jacobian():
    hill, grads = _stub_hill()
    ud = _StubUD(orography=hill, orography_grad=grads)
    elem, _ = dis_grid.grid_init(ud)

    class _Folded(_Tier2StretchMap):
        def tangents(self, xi):
            t = super().tangents(xi)
            t[0][0] = t[0][0] - 2.0  # x' < 0: orientation flips
            return t

    bad = _Folded(hill, grads, eta0=ud.ymin, etat=ud.ymax)
    with pytest.raises(ValueError, match="Jacobian"):
        terrain.build_metric_fields_from_map(elem, ud, bad)
