"""Phase-0 reduction contract for the curvilinear metric generalization.

Executable definition of "reduces to the current terrain code" for the
planned Klein-style metric (J, N1, N2, N3): for the vertical-line map
z = z(xi_h1, eta, xi_h2) (x, y identity) the tangents in role order
(h1, v, h2) are

    t1 = (1, G1, 0),   t2 = (0, J, 0),   t3 = (0, G2, 1)

with J = dz/deta and G_h = dz/dxi_h, and the area normals N_i = t_j x t_k
(cyclic) must come out as

    N1 = (J, 0, 0),   N2 = (-G1, 1, -G2),   N3 = (0, 0, J)

so the general flux components F_i = N_i . f reproduce today's
``_metric_contravariant_fluxes_jit`` BIT-EXACTLY:

    F1 = J f_h1,   F2 = f_v - G1 f_h1 - G2 f_h2,   F3 = J f_h2.

The bit-exact contraction order is part of the contract: vertical
component first, then h1, then h2 (see ``general_flux``) — Phase 2's
general flux assembly must keep this order so the reduction stays exact.

Also asserted: the duality N_i . t_j = J delta_ij (equivalent to
N_i = J grad xi_i, the bridge to the gradient/elliptic picture),
det N = J^2 (nonsingularity), and a Tier-2 fixture
x = (x(xi1), z(xi1, eta, xi3), y(xi3)) pinning the identities Phase 4's
stretched-grid gate depends on — in particular the effective-slope
identity -(N2)_h / (N2)_v = z_xi_h / x'_h that keeps the contravariant
wall reflection (cell_boundary._slope_terms) correct on stretched grids.
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
    """The Phase-2 contraction contract: F_i = N_i . f, vertical first.

    Role-component order is (h1, v, h2); summing v + h1 + h2 left-to-right
    makes the vertical-line reduction bit-exact against the legacy kernel
    (1*f_v == f_v, x + (-y) == x - y, 0 + x == x in IEEE arithmetic).
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


# ----------------------------------------- vertical-line maps (Tier 1)


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


# ------------------------------------------------ Tier-2 fixture (stretch)


def _tier2_map():
    """Analytic Tier-2 map x = (x(xi1), z(xi1, eta, xi3), y(xi3)).

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

    # closed forms of the Tier-2 normals and Jacobian
    np.testing.assert_allclose(J, xp * z_eta * yp, rtol=1e-14)
    assert np.all(J > 0.0)
    np.testing.assert_allclose(N1[0], z_eta * yp, rtol=1e-14)
    assert np.array_equal(N1[1], np.zeros_like(J))
    np.testing.assert_allclose(N2[0], -(yp * z_x1), rtol=1e-14)
    np.testing.assert_allclose(N2[1], xp * yp, rtol=1e-14)
    np.testing.assert_allclose(N2[2], -(xp * z_x3), rtol=1e-14)
    np.testing.assert_allclose(N3[2], xp * z_eta, rtol=1e-14)

    # duality and nonsingularity at Tier 2 (now with real FP cancellation)
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
