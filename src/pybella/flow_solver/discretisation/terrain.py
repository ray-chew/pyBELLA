"""Terrain-following vertical coordinates — transform + metric fields.

The solver integrates in computational coordinates ``(xi_h1, eta, xi_h2)``
(role order: first horizontal, vertical, second horizontal) on the existing
uniform grid; terrain enters only through precomputed metric fields. With
terrain height ``h(xi_h1, xi_h2)`` and the flat domain top at ``eta = etat``,
a :class:`VerticalTransform` defines the physical height ``z(xi, eta)`` and
the two metric quantities every operator consumes:

    J   = dz/deta                  (Jacobian; cell "thickness" weight)
    G_i = dz/dxi_i at fixed eta    (slope terms, i in {h1, h2})

Gal-Chen--Somerville (:class:`GalChenTransform`) is the first concrete
transform; SLEVE drops in later as another subclass — all metric arrays are
stored as full grid-shaped fields even where Gal-Chen makes them separable,
so no operator changes are needed for an eta-dependent Jacobian.

Activation contract (h == 0 bypass): terrain is active iff ``ud.orography``
is defined (a callable ``h(xi_h1, xi_h2)`` in nondimensional units, mirroring
the ``ud.stratification`` convention). Without it :func:`build_metric_fields`
returns ``None``, ``elem.metric``/``node.metric`` are ``None``, and every
consumer takes the uniform-Cartesian code path untouched — bit-identity with
the pre-terrain solver by construction.

Slopes prefer an analytic gradient ``ud.orography_grad = (dh_dxi1, dh_dxi2)``
(role-ordered callables); otherwise they are central differences of the
orography callable evaluated at shifted coordinates (exact at ghost cells,
no array stenciling).

Everything is role-oriented through :mod:`pybella.utils.axes` — the metric
machinery works for any ``ud.gravity_direction`` in {0, 1, 2}.
"""

import numpy as np

from ...utils import axes
from ...utils import options as opts


class VerticalTransform:
    """Abstract map from computational vertical eta in [eta0, etat] to z.

    Subclasses implement elementwise, broadcastable methods; the builder
    materialises full grid-shaped arrays from them.
    """

    def z(self, eta, h, eta0, etat):
        """Physical height z(eta, h)."""
        raise NotImplementedError

    def jacobian(self, eta, h, eta0, etat):
        """dz/deta at fixed horizontal position."""
        raise NotImplementedError

    def decay(self, eta, eta0, etat):
        """Weight b(eta) by which terrain influence decays with height."""
        raise NotImplementedError

    def slope(self, eta, dh, eta0, etat):
        """G = dz/dxi at fixed eta, given the terrain slope dh = dh/dxi.

        Default assumes ``z = eta + h * b(eta)``; transforms with a
        different structure (e.g. SLEVE's two-scale split) override this.
        """
        return dh * self.decay(eta, eta0, etat)


class GalChenTransform(VerticalTransform):
    """Gal-Chen & Somerville (1975): linear decay of terrain with height.

    z = eta + h * b(eta),  b(eta) = (etat - eta) / (etat - eta0)

    so z(eta0) = eta0 + h (surface follows the terrain) and z(etat) = etat
    (flat top). The Jacobian J = 1 - h / (etat - eta0) is eta-independent.
    """

    def decay(self, eta, eta0, etat):
        return (etat - eta) / (etat - eta0)

    def z(self, eta, h, eta0, etat):
        return eta + h * self.decay(eta, eta0, etat)

    def jacobian(self, eta, h, eta0, etat):
        # broadcast against eta so the builder always gets a full field
        return (1.0 - h / (etat - eta0)) + 0.0 * eta


class MetricFields:
    """Precomputed terrain metric arrays on one grid (cells or nodes).

    Arrays are array-axis oriented (same layout as the solution fields,
    squeezed like :class:`CellSolField`); the slope arrays ``G1``/``G2``
    are the dz/dxi terms along the horizontal role axes ``haxes`` mapped
    by ``axes.role_perm``. ``G2`` is ``None`` in 2D. Plain float64 arrays
    only, safe to pass straight into numba kernels.
    """

    _ARRAYS = ("J", "ooJ", "G1", "G2", "z")

    def __init__(self, J, G1, G2, z, vaxis, haxes):
        self.J = J
        self.ooJ = 1.0 / J
        self.G1 = G1
        self.G2 = G2
        self.z = z
        self.vaxis = vaxis
        self.haxes = haxes

    def flip_forward(self):
        """Mirror CellSolField.flip_forward for the advection sweeps."""
        for key in self._ARRAYS:
            value = getattr(self, key)
            if value is not None:
                setattr(self, key, np.moveaxis(value, 0, -1))
        self.vaxis, self.haxes = self._shift_axes(-1)

    def flip_backward(self):
        for key in self._ARRAYS:
            value = getattr(self, key)
            if value is not None:
                setattr(self, key, np.moveaxis(value, -1, 0))
        self.vaxis, self.haxes = self._shift_axes(+1)

    def _shift_axes(self, step):
        ndim = self.J.ndim
        shift = lambda a: (a + step) % ndim if a is not None else None
        return shift(self.vaxis), tuple(shift(a) for a in self.haxes)


def apply_gradient_map(metric, dp):
    """Physical gradients from computational ones: dp <- A @ dp (in place).

    ``dp`` is the axis-indexed list of the three cell-gradient arrays. The
    chain rule for z = z(xi, eta) gives, in role order (h1, v, h2),

        d/dx_h|z = d/dxi_h - (G_h / J) d/deta,   d/dz = (1/J) d/deta

    i.e. the matrix A = [[1, -G1/J, 0], [0, 1/J, 0], [0, -G2/J, 1]]. The
    horizontal rows are corrected before the vertical row is scaled.
    """
    a_h1, a_h2 = metric.haxes
    dp_v = dp[metric.vaxis]
    dp[a_h1] = dp[a_h1] - metric.G1 * metric.ooJ * dp_v
    if a_h2 is not None:
        dp[a_h2] = dp[a_h2] - metric.G2 * metric.ooJ * dp_v
    dp[metric.vaxis] = dp_v * metric.ooJ
    return dp


def elliptic_tensor(metric, h_role):
    """Fold the terrain metric into the role-indexed H^-1 tensor.

    Returns M = J A^T H^-1 A (role order h1, v, h2), the coefficient
    tensor of the elliptic operator: the rhs divergence measures
    D_i((J A^T F)_i) and the momentum correction applies H^-1 A grad p,
    so their composition carries exactly this tensor. With H^-1 == I it
    is the classic terrain-following tensor

        [[J, -G1, 0], [-G1, (1 + G1^2 + G2^2)/J, -G2], [0, -G2, J]]

    and with h == 0 (J == 1, G == 0) it reduces bit-exactly to ``h_role``.
    """
    J, ooJ, G1, G2 = metric.J, metric.ooJ, metric.G1, metric.G2
    h = h_role
    M00 = J * h[0][0]
    M02 = J * h[0][2]
    M20 = J * h[2][0]
    M22 = J * h[2][2]
    M01 = -G1 * h[0][0] + h[0][1] - G2 * h[0][2]
    M10 = -G1 * h[0][0] + h[1][0] - G2 * h[2][0]
    M12 = -G1 * h[0][2] + h[1][2] - G2 * h[2][2]
    M21 = -G1 * h[2][0] + h[2][1] - G2 * h[2][2]
    M11 = ooJ * (
        G1 * G1 * h[0][0]
        - G1 * (h[0][1] + h[1][0])
        + G1 * G2 * (h[0][2] + h[2][0])
        + h[1][1]
        - G2 * (h[1][2] + h[2][1])
        + G2 * G2 * h[2][2]
    )
    return ((M00, M01, M02), (M10, M11, M12), (M20, M21, M22))


def terrain_is_active(ud):
    return getattr(ud, "orography", None) is not None


def get_transform(ud):
    transform = getattr(ud, "vertical_transform", None)
    return transform if transform is not None else GalChenTransform()


def vertical_extent(ud, v):
    """(eta0, etat): domain extent along the vertical axis v."""
    return ((ud.xmin, ud.xmax), (ud.ymin, ud.ymax), (ud.zmin, ud.zmax))[v]


def _coordinate_wrap(ud, axis):
    """Identity, or periodic wrap into the domain for PERIODIC axes.

    Ghost coordinates lie outside the domain; every other field sees its
    periodic image there (ghost-cell wrap), so the orography must too —
    otherwise the metric is discontinuous across the periodic seam and the
    elliptic system becomes inconsistent at the duplicated nodes.
    """
    if axis is None or ud.bdry_type[axis] != opts.BdryType.PERIODIC:
        return lambda c: c
    lo, hi = vertical_extent(ud, axis)
    length = hi - lo
    return lambda c: lo + np.mod(c - lo, length)


def _effective_orography(ud, a_h1, a_h2):
    """ud.orography with periodic-wrapped arguments (h1, h2 role order)."""
    wrap1 = _coordinate_wrap(ud, a_h1)
    wrap2 = _coordinate_wrap(ud, a_h2)
    return lambda xi1, xi2: ud.orography(wrap1(xi1), wrap2(xi2))


def _coord_view(grid_obj, axis, ndim):
    """Coordinate array of `axis`, shaped to broadcast over an ndim field."""
    shape = [1] * ndim
    shape[axis] = -1
    return axes.coords_along(grid_obj, axis).reshape(shape)


def _terrain_slope(ud, heff, a_h1, a_h2, xi1, xi2, which, spacing):
    """dh/dxi_which (role index 0 or 1): analytic if provided, else FD.

    Both paths wrap periodic coordinates: the analytic gradient is
    evaluated at the wrapped points, the central difference differentiates
    the wrapped (periodic) effective orography so the seam is consistent.
    """
    grad = getattr(ud, "orography_grad", None)
    if grad is not None:
        wrap1 = _coordinate_wrap(ud, a_h1)
        wrap2 = _coordinate_wrap(ud, a_h2)
        return grad[which](wrap1(xi1), wrap2(xi2))
    d = spacing
    if which == 0:
        return (heff(xi1 + d, xi2) - heff(xi1 - d, xi2)) / (2.0 * d)
    return (heff(xi1, xi2 + d) - heff(xi1, xi2 - d)) / (2.0 * d)


def build_metric_fields(grid_obj, ud):
    """Build MetricFields for one grid (ElemSpaceDiscr or NodeSpaceDiscr).

    Returns None when terrain is inactive — callers branch on that and the
    uniform-Cartesian path stays untouched.
    """
    if not terrain_is_active(ud):
        return None

    ndim = grid_obj.ndim
    v = axes.vertical_axis(ud)
    transform = get_transform(ud)
    eta0, etat = vertical_extent(ud, v)

    if ndim == 2:
        # axes.validate enforces v == 1 in 2D: x horizontal, no second
        # horizontal (terrain runs are quasi-2D 3D for now, but the metric
        # build supports native 2D for the planned lap2D cross-term work)
        a_h1, a_h2 = 0, None
    else:
        a_h1, a_h2 = axes.horizontal_axes(v)

    shape = tuple(int(grid_obj.sc[dim]) for dim in range(ndim))
    eta = _coord_view(grid_obj, v, ndim)
    xi1 = _coord_view(grid_obj, a_h1, ndim)
    xi2 = _coord_view(grid_obj, a_h2, ndim) if a_h2 is not None else 0.0

    heff = _effective_orography(ud, a_h1, a_h2)
    h = heff(xi1, xi2)

    def full(expr):
        return np.ascontiguousarray(
            np.broadcast_to(expr, shape).astype(np.float64, copy=False)
        )

    J = full(transform.jacobian(eta, h, eta0, etat))
    if np.any(J <= 0.0):
        raise ValueError(
            "terrain transform produced non-positive Jacobian: "
            "orography reaches or exceeds the domain top"
        )

    z = full(transform.z(eta, h, eta0, etat))

    dh1 = _terrain_slope(ud, heff, a_h1, a_h2, xi1, xi2, 0, grid_obj.dxyz[a_h1])
    G1 = full(transform.slope(eta, dh1, eta0, etat))
    if a_h2 is not None:
        dh2 = _terrain_slope(ud, heff, a_h1, a_h2, xi1, xi2, 1, grid_obj.dxyz[a_h2])
        G2 = full(transform.slope(eta, dh2, eta0, etat))
    else:
        G2 = None

    return MetricFields(J=J, G1=G1, G2=G2, z=z, vaxis=v, haxes=(a_h1, a_h2))
