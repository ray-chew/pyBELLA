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

    def __init__(self, v=1, orography=None, orography_grad=None):
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
