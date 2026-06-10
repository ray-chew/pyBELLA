"""Unit tests for the terrain-following transform + metric field builder.

Covers the Gal-Chen analytic identities, the SLEVE-ready interface
contract (slope = decay * dh), the h == 0 exact-identity guarantee, and
the role orientation of the metric arrays for every gravity_direction.
"""

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import terrain
from pybella.utils import axes
from pybella.utils import options as opts

# --- transform-level identities --------------------------------------------


def test_galchen_surface_and_top():
    tr = terrain.GalChenTransform()
    eta0, etat = 0.0, 2.0
    h = np.array([0.0, 0.1, 0.3])
    # surface follows terrain, top is flat
    assert np.allclose(tr.z(eta0, h, eta0, etat), eta0 + h)
    assert np.allclose(tr.z(etat, h, eta0, etat), etat)


def test_galchen_jacobian_matches_fd():
    tr = terrain.GalChenTransform()
    eta0, etat = 0.5, 3.0
    h = 0.2
    eta = np.linspace(eta0, etat, 11)
    d = 1e-6
    fd = (tr.z(eta + d, h, eta0, etat) - tr.z(eta - d, h, eta0, etat)) / (2 * d)
    assert np.allclose(tr.jacobian(eta, h, eta0, etat), fd, atol=1e-9)
    # Gal-Chen J is eta-independent and equals 1 - h/(etat - eta0)
    assert np.allclose(tr.jacobian(eta, h, eta0, etat), 1.0 - h / (etat - eta0))


def test_galchen_slope_is_decay_weighted():
    tr = terrain.GalChenTransform()
    eta0, etat = 0.0, 1.0
    eta = np.linspace(eta0, etat, 5)
    dh = 0.07
    assert np.allclose(tr.slope(eta, dh, eta0, etat), dh * (etat - eta) / (etat - eta0))
    # slope vanishes at the flat top, equals dh at the surface
    assert tr.slope(etat, dh, eta0, etat) == 0.0
    assert tr.slope(eta0, dh, eta0, etat) == dh


# --- builder ----------------------------------------------------------------


class _StubUD:
    """Minimal ud for grid_init + build_metric_fields."""

    def __init__(
        self,
        v=1,
        orography=None,
        orography_grad=None,
        orography_smooth=None,
        orography_smooth_grad=None,
        vertical_transform=None,
    ):
        self.inx, self.iny, self.inz = 9, 7, 5
        self.xmin, self.xmax = -1.0, 1.0
        self.ymin, self.ymax = 0.0, 2.0
        self.zmin, self.zmax = -0.5, 0.5
        self.bdry_type = np.array(
            [opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.PERIODIC]
        )
        self.gravity_direction = v
        if orography is not None:
            self.orography = orography
        if orography_grad is not None:
            self.orography_grad = orography_grad
        if orography_smooth is not None:
            self.orography_smooth = orography_smooth
        if orography_smooth_grad is not None:
            self.orography_smooth_grad = orography_smooth_grad
        if vertical_transform is not None:
            self.vertical_transform = vertical_transform


def _hill(amplitude=0.05, width=0.3):
    def h(xi1, xi2):
        return amplitude / (1.0 + (xi1 / width) ** 2) + 0.0 * xi2

    def dh1(xi1, xi2):
        return (
            -2.0 * amplitude * xi1 / width**2 / (1.0 + (xi1 / width) ** 2) ** 2
            + 0.0 * xi2
        )

    def dh2(xi1, xi2):
        return 0.0 * xi1 + 0.0 * xi2

    return h, (dh1, dh2)


def test_metric_none_without_orography():
    elem, node = dis_grid.grid_init(_StubUD())
    assert elem.metric is None
    assert node.metric is None


def test_h_zero_gives_exact_identity_metric():
    flat = lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2
    ud = _StubUD(orography=flat)
    elem, node = dis_grid.grid_init(ud)
    for grid_obj in (elem, node):
        m = grid_obj.metric
        assert m is not None
        v = axes.vertical_axis(ud)
        # exact identities, not approximate: h == 0 incurs no rounding
        assert np.all(m.J == 1.0)
        assert np.all(m.ooJ == 1.0)
        assert np.all(m.G1 == 0.0)
        assert np.all(m.G2 == 0.0)
        eta = axes.coords_along(grid_obj, v)
        shape = [1] * grid_obj.ndim
        shape[v] = -1
        assert np.all(m.z == np.broadcast_to(eta.reshape(shape), m.z.shape))


@pytest.mark.parametrize("v", [0, 1, 2])
def test_role_orientation(v):
    h, grad = _hill()
    ud = _StubUD(v=v, orography=h, orography_grad=grad)
    elem, _ = dis_grid.grid_init(ud)
    m = elem.metric
    a_h1, a_h2 = axes.horizontal_axes(v)
    assert m.vaxis == v
    assert m.haxes == (a_h1, a_h2)
    assert m.J.shape == tuple(int(elem.sc[d]) for d in range(elem.ndim))

    # J varies along h1 (the hill axis) and is constant along v and h2
    assert np.ptp(m.J, axis=a_h1).max() > 0.0
    assert np.ptp(m.J, axis=v).max() == 0.0  # Gal-Chen: eta-independent
    assert np.ptp(m.J, axis=a_h2).max() == 0.0

    # G1 decays with height: max slope at the bottom layer, ~0 at the top
    lo, hi = axes.wall_slabs(elem.ndim, v, depth=1)
    assert np.abs(m.G1[lo]).max() > np.abs(m.G1[hi]).max()
    # z is strictly increasing along the vertical
    assert np.all(np.diff(m.z, axis=v) > 0.0)
    # flat ridge: no slope along h2
    assert np.all(m.G2 == 0.0)


def test_fd_slope_matches_analytic_grad():
    h, grad = _hill()
    ud_fd = _StubUD(orography=h)
    ud_an = _StubUD(orography=h, orography_grad=grad)
    # resolve the hill (width 0.3) properly so the second-order FD converges;
    # walls in x: the hill does not decay to zero at the domain edge, so the
    # periodic coordinate wrap would (correctly) introduce a seam kink that
    # FD smears but the analytic gradient does not — not what's tested here
    for ud in (ud_fd, ud_an):
        ud.inx = 129
        ud.bdry_type = np.array([opts.BdryType.WALL] * 3)
    elem_fd, _ = dis_grid.grid_init(ud_fd)
    elem_an, _ = dis_grid.grid_init(ud_an)
    # ~4e-4 truncation at the steepest point on this grid; the test guards
    # orientation/sign/scale, not FD order
    assert np.allclose(elem_fd.metric.G1, elem_an.metric.G1, atol=1e-3)


def test_orography_above_top_rejected():
    tall = lambda xi1, xi2: 2.5 + 0.0 * xi1 + 0.0 * xi2  # > ymax - ymin
    with pytest.raises(ValueError, match="Jacobian"):
        dis_grid.grid_init(_StubUD(orography=tall))


def test_flip_forward_backward_roundtrip():
    h, grad = _hill()
    ud = _StubUD(orography=h, orography_grad=grad)
    elem, _ = dis_grid.grid_init(ud)
    m = elem.metric
    J0, G10, vax0, hax0 = m.J.copy(), m.G1.copy(), m.vaxis, m.haxes
    m.flip_forward()
    assert m.J.shape == tuple(np.roll(J0.shape, -1))
    assert m.vaxis == (vax0 - 1) % elem.ndim
    m.flip_backward()
    assert np.array_equal(m.J, J0)
    assert np.array_equal(m.G1, G10)
    assert m.vaxis == vax0
    assert m.haxes == hax0


# --- SLEVE -------------------------------------------------------------------


def _sleve(n=1.0):
    # decay scales well separated and inside the [0, 2] stub domain
    return terrain.SLEVETransform(s1=1.5, s2=0.25, n=n)


def _ridge(amplitude=0.04, k=12.0):
    """Small-scale ridge h2 on a smooth envelope h1 (analytic split)."""

    def h1(xi1, xi2):
        return amplitude * np.exp(-(xi1**2)) + 0.0 * xi2

    def h(xi1, xi2):
        return h1(xi1, xi2) * (1.0 + 0.5 * np.cos(k * xi1))

    def dh1_1(xi1, xi2):
        return -2.0 * xi1 * amplitude * np.exp(-(xi1**2)) + 0.0 * xi2

    def dh_1(xi1, xi2):
        return dh1_1(xi1, xi2) * (1.0 + 0.5 * np.cos(k * xi1)) - h1(
            xi1, xi2
        ) * 0.5 * k * np.sin(k * xi1)

    zero = lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2
    return h, (dh_1, zero), h1, (dh1_1, zero)


@pytest.mark.parametrize("n", [1.0, 1.35])
def test_sleve_surface_and_top(n):
    tr = _sleve(n)
    eta0, etat = 0.0, 2.0
    h1 = np.array([0.0, 0.05, 0.1])
    h2 = np.array([0.0, 0.03, -0.02])
    z_surf = tr.z(eta0, (h1, h2), eta0, etat)
    assert np.allclose(z_surf, eta0 + h1 + h2)
    assert np.allclose(tr.z(etat, (h1, h2), eta0, etat), etat)


@pytest.mark.parametrize("n", [1.0, 1.35])
def test_sleve_jacobian_matches_fd(n):
    tr = _sleve(n)
    eta0, etat = 0.0, 2.0
    h = (0.08, 0.04)
    # stay inside (eta0, etat): the FD stencil must not straddle the
    # surface where the below-ground linear extension kicks in
    eta = np.linspace(eta0 + 1e-3, etat - 1e-3, 41)
    d = 1e-7
    fd = (tr.z(eta + d, h, eta0, etat) - tr.z(eta - d, h, eta0, etat)) / (2 * d)
    assert np.allclose(tr.jacobian(eta, h, eta0, etat), fd, atol=1e-6)


def test_sleve_jacobian_eta_dependent():
    tr = _sleve()
    eta0, etat = 0.0, 2.0
    eta = np.linspace(eta0, etat, 9)
    J = tr.jacobian(eta, (0.08, 0.04), eta0, etat)
    assert np.ptp(J) > 1e-3  # the first eta-dependent Jacobian


def test_sleve_n_gt_1_has_uniform_surface_jacobian():
    tr = _sleve(n=1.35)
    eta0, etat = 0.0, 2.0
    J_surf = tr.jacobian(eta0, (0.08, 0.04), eta0, etat)
    assert np.allclose(J_surf, 1.0)  # db_i(0) = 0 for n > 1


@pytest.mark.parametrize("n", [1.0, 1.35])
def test_sleve_below_ground_extension_is_finite_and_smooth(n):
    tr = _sleve(n)
    eta0, etat = 0.0, 2.0
    eta_ghost = np.array([-0.4, -0.2, -1e-9])
    for fn in (tr.z, tr.jacobian):
        vals = fn(eta_ghost, (0.08, 0.04), eta0, etat)
        assert np.all(np.isfinite(vals))
    # continuity across the surface (for n > 1 the decay slope behaves as
    # zeta^(n-1) just above ground — continuous but steep, hence the loose
    # tolerance; n == 1 matches to machine precision)
    assert np.allclose(
        tr.jacobian(-1e-9, (0.08, 0.04), eta0, etat),
        tr.jacobian(+1e-9, (0.08, 0.04), eta0, etat),
        atol=1e-3 if n > 1.0 else 1e-9,
    )


def test_sleve_scale_separation():
    """The point of SLEVE: at mid-levels the small-scale decay b2 is far
    below b1, so pure small-scale terrain barely distorts the grid there —
    under Gal-Chen it still carries ~half its surface slope."""
    tr = _sleve()
    gc = terrain.GalChenTransform()
    eta0, etat = 0.0, 2.0
    eta_mid = 1.0
    dh = 1.0  # pure small-scale slope
    g_sleve = tr.slope(eta_mid, (0.0, dh), eta0, etat)
    g_gc = gc.slope(eta_mid, dh, eta0, etat)
    assert abs(g_sleve) < 0.05 * abs(g_gc)
    # and b2 <= b1 everywhere for s2 < s1
    eta = np.linspace(eta0, etat, 33)
    b1, _ = tr._b_db(eta, eta0, etat, tr.s1)
    b2, _ = tr._b_db(eta, eta0, etat, tr.s2)
    assert np.all(b2 <= b1 + 1e-12)


def test_sleve_builder_requires_smooth_split():
    h, grad, _, _ = _ridge()
    ud = _StubUD(orography=h, orography_grad=grad, vertical_transform=_sleve())
    with pytest.raises(ValueError, match="orography_smooth"):
        dis_grid.grid_init(ud)


def test_sleve_h_zero_gives_exact_identity_metric():
    flat = lambda xi1, xi2: 0.0 * xi1 + 0.0 * xi2
    ud = _StubUD(orography=flat, orography_smooth=flat, vertical_transform=_sleve())
    elem, node = dis_grid.grid_init(ud)
    for grid_obj in (elem, node):
        m = grid_obj.metric
        assert np.all(m.J == 1.0)
        assert np.all(m.G1 == 0.0)
        assert np.all(m.G2 == 0.0)


def test_sleve_builder_slope_matches_fd_of_z():
    """G1 from the builder == d z / d xi1 at fixed eta, by FD across columns."""
    h, grad, h1, grad1 = _ridge()
    ud = _StubUD(
        orography=h,
        orography_grad=grad,
        orography_smooth=h1,
        orography_smooth_grad=grad1,
        vertical_transform=_sleve(),
    )
    ud.inx = 257  # resolve the k = 12 ridge for the cross-column FD
    elem, _ = dis_grid.grid_init(ud)
    m = elem.metric
    dz_dxi1 = np.gradient(m.z, elem.dx, axis=0)
    inner = (slice(4, -4), slice(None), slice(None))
    assert np.allclose(m.G1[inner], dz_dxi1[inner], atol=2e-3)


def test_sleve_jacobian_positivity_rejected():
    # small-scale amplitude ~ s2: the residual decay overshoots J <= 0
    h, grad, h1, grad1 = _ridge(amplitude=0.3, k=12.0)
    ud = _StubUD(
        orography=h,
        orography_grad=grad,
        orography_smooth=h1,
        orography_smooth_grad=grad1,
        vertical_transform=terrain.SLEVETransform(s1=1.5, s2=0.05),
    )
    with pytest.raises(ValueError, match="Jacobian"):
        dis_grid.grid_init(ud)


def test_sleve_flip_roundtrip_with_eta_dependent_J():
    h, grad, h1, grad1 = _ridge()
    ud = _StubUD(
        orography=h,
        orography_grad=grad,
        orography_smooth=h1,
        orography_smooth_grad=grad1,
        vertical_transform=_sleve(),
    )
    elem, _ = dis_grid.grid_init(ud)
    m = elem.metric
    assert np.ptp(m.J, axis=m.vaxis).max() > 0.0  # genuinely eta-dependent
    J0, G10, vax0, hax0 = m.J.copy(), m.G1.copy(), m.vaxis, m.haxes
    m.flip_forward()
    m.flip_backward()
    assert np.array_equal(m.J, J0)
    assert np.array_equal(m.G1, G10)
    assert (m.vaxis, m.haxes) == (vax0, hax0)
