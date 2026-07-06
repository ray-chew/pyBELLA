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

    ``n_components`` declares how many orography fields the transform
    consumes. Single-component transforms (the default) receive plain
    arrays for ``h``/``dh``; two-component transforms (SLEVE) receive
    tuples ``h = (h_smooth, h_residual)`` / ``dh = (dh_smooth,
    dh_residual)`` — the split is built from ``ud.orography_smooth``
    against the total ``ud.orography`` by :func:`build_metric_fields`.
    """

    n_components = 1

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


class SLEVETransform(VerticalTransform):
    """SLEVE (Schär et al. 2002; Leuenberger et al. 2010 exponent n).

    Two-scale split z = eta + h1 b1(eta) + h2 b2(eta) with per-component
    decay (zeta = eta - eta0, H = etat - eta0)

        b_i(zeta) = sinh((H/s_i)^n - (zeta/s_i)^n) / sinh((H/s_i)^n)

    so b_i(0) = 1 (terrain-following surface), b_i(H) = 0 (flat top) and
    the small-scale part h2 decays on its own scale s2 << s1 — the point
    of SLEVE: small-scale terrain distortion leaves the grid quickly with
    height instead of propagating to every level as under Gal-Chen.

    ``s1``/``s2`` are nondimensional decay heights (same units as eta);
    ``n > 1`` (e.g. Leuenberger's 1.35) gives db_i(0) = 0, i.e. an exactly
    uniform Jacobian at the surface — the steep-terrain fallback. Below
    the surface (ghost rows, zeta < 0) the decay continues linearly with
    its surface slope, which keeps fractional ``n`` well-defined and J
    smooth across the bottom boundary.

    The first transform with an eta-dependent Jacobian — all operators
    consume full-field J(eta) so nothing downstream changes.
    """

    n_components = 2

    def __init__(self, s1, s2, n=1.0):
        self.s1 = float(s1)
        self.s2 = float(s2)
        self.n = float(n)

    def _b_db(self, eta, eta0, etat, s):
        """Decay b and its eta-derivative db, linearly extended below eta0."""
        n = self.n
        H = etat - eta0
        zeta = eta - eta0 + 0.0 * np.asarray(eta)
        zc = np.maximum(zeta, 0.0)
        arg_top = (H / s) ** n
        oosinh = 1.0 / np.sinh(arg_top)
        inner = arg_top - (zc / s) ** n
        b = np.sinh(inner) * oosinh
        db = -(n * zc ** (n - 1.0) / s**n) * np.cosh(inner) * oosinh
        if n == 1.0:
            db0 = -(1.0 / s) * np.cosh(arg_top) * oosinh
        else:
            db0 = 0.0  # n > 1: zero surface slope of the decay
        b = np.where(zeta < 0.0, 1.0 + db0 * zeta, b)
        db = np.where(zeta < 0.0, db0 + 0.0 * zeta, db)
        return b, db

    def z(self, eta, h, eta0, etat):
        h1, h2 = h
        b1, _ = self._b_db(eta, eta0, etat, self.s1)
        b2, _ = self._b_db(eta, eta0, etat, self.s2)
        return eta + h1 * b1 + h2 * b2

    def jacobian(self, eta, h, eta0, etat):
        h1, h2 = h
        _, db1 = self._b_db(eta, eta0, etat, self.s1)
        _, db2 = self._b_db(eta, eta0, etat, self.s2)
        return 1.0 + h1 * db1 + h2 * db2

    def slope(self, eta, dh, eta0, etat):
        # two decays — the single-decay base default cannot express this
        dh1, dh2 = dh
        b1, _ = self._b_db(eta, eta0, etat, self.s1)
        b2, _ = self._b_db(eta, eta0, etat, self.s2)
        return dh1 * b1 + dh2 * b2


class MetricFields:
    """Precomputed terrain metric arrays on one grid (cells or nodes).

    Arrays are array-axis oriented (same layout as the solution fields,
    squeezed like :class:`CellSolField`); the slope arrays ``G1``/``G2``
    are the dz/dxi terms along the horizontal role axes ``haxes`` mapped
    by ``axes.role_perm``. ``G2`` is ``None`` in 2D. Plain float64 arrays
    only, safe to pass straight into numba kernels.

    General (Klein) metric data rides alongside the legacy scalars:

    ``N``
        Area normals N_a = t_b x t_c (cyclic over computational axes).
        ``N[a]`` is the normal of the xi_a = const surface for the CURRENT
        array axis ``a`` (the outer index rotates with the sweep flips);
        its entries are the ndim CARTESIAN components as grid-shaped
        arrays (the inner index is fixed — momenta keep their identity
        under flips, so component k always multiplies momentum k).
    ``x``
        Physical coordinates, Cartesian-indexed; ``x[cart_v]`` is the
        height field (generalizes ``z``), entries are ``None`` where the
        map is the identity along that axis (not materialized).
    ``cart_v`` / ``cart_haxes``
        The fixed Cartesian role axes (vaxis/haxes at build time); unlike
        ``vaxis``/``haxes`` they never change under flips.

    When ``N``/``x`` are not supplied they are synthesized from the
    vertical-line map's normals (J, 0, 0), (-G1, 1, -G2), (0, 0, J) —
    bit-identical to the cross products of the tangents
    t1 = (1, G1, 0), t2 = (0, J, 0), t3 = (0, G2, 1) (the Phase-0
    reduction contract, ``test_scripts/test_metric_reduction.py``).

    Beyond-vertical-line (Tier 3) metric data, all inert on the legacy
    path:

    ``height``
        Generalized altitude (the coordinate gravity acts along). For
        vertical-line maps this IS ``z`` (aliased); a curved map (sphere)
        supplies e.g. ``r - a``. Consumers that mean "height above the
        reference geopotential" read this, not ``z``.
    ``h_v``
        Vertical arc length per unit eta, |t_v|. For vertical-line maps
        this IS ``J`` (aliased; t_v = (0, J, 0)), and it equals the
        legacy ghost-spacing construction J / (N_v)_v there.
    ``e_up``
        Unit "up" direction as a Cartesian-component list (like one row
        of ``N``), or ``None`` when up is the fixed Cartesian role-v axis
        (every vertical-line map). Buoyancy/gravity consumers branch on
        this.
    ``vertical_line``
        False for maps whose vertical coordinate lines are not parallel
        Cartesian lines. Then ``G1``/``G2`` are ``None`` — the slope
        scalars are mathematically undefined (their construction divides
        by a Cartesian component of N_v that passes through zero on a
        sphere) — and any remaining G1/G2 consumer must not be reached.
    """

    _ARRAYS = ("J", "ooJ", "G1", "G2", "z", "height", "h_v")

    def __init__(
        self,
        J,
        G1,
        G2,
        z,
        vaxis,
        haxes,
        N=None,
        x=None,
        height=None,
        h_v=None,
        e_up=None,
        vertical_line=True,
    ):
        self.J = J
        self.ooJ = 1.0 / J
        self.G1 = G1
        self.G2 = G2
        self.z = z
        self.vaxis = vaxis
        self.haxes = haxes
        self.cart_v = vaxis
        self.cart_haxes = haxes
        self.N = self._vertical_line_normals() if N is None else N
        if x is None:
            x = [None] * J.ndim
            x[vaxis] = z
        self.x = x
        self.height = z if height is None else height
        self.h_v = J if h_v is None else h_v
        self.e_up = e_up
        self.vertical_line = vertical_line

    def _vertical_line_normals(self):
        """Cartesian-component normals of the vertical-line map (canonical
        orientation: array axes == Cartesian axes at construction)."""
        ndim = self.J.ndim
        one = np.ones_like(self.J)
        zero = np.zeros_like(self.J)
        a_h1, a_h2 = self.haxes
        v = self.vaxis
        N = [[zero] * ndim for _ in range(ndim)]
        # N_h1 = (J along cart h1); N_v = (-G_h horizontals, 1 vertical)
        N[a_h1][a_h1] = self.J
        N[v][a_h1] = -self.G1
        N[v][v] = one
        if a_h2 is not None:
            N[v][a_h2] = -self.G2
            N[a_h2][a_h2] = self.J
        return N

    def flip_forward(self):
        """Mirror CellSolField.flip_forward for the advection sweeps.

        Every leaf array rolls its axes; the OUTER index of ``N`` rotates
        with them (the normal of old array axis a lands on (a - 1) % ndim)
        while the Cartesian component index never moves. ``cart_*`` are
        flip-invariant by definition.
        """
        roll = lambda arr: np.moveaxis(arr, 0, -1)
        for key in self._ARRAYS:
            value = getattr(self, key)
            if value is not None:
                setattr(self, key, roll(value))
        rolled = [[roll(c) for c in Na] for Na in self.N]
        self.N = rolled[1:] + rolled[:1]
        self.x = [None if c is None else roll(c) for c in self.x]
        if self.e_up is not None:
            self.e_up = [roll(c) for c in self.e_up]
        self.vaxis, self.haxes = self._shift_axes(-1)

    def flip_backward(self):
        roll = lambda arr: np.moveaxis(arr, -1, 0)
        for key in self._ARRAYS:
            value = getattr(self, key)
            if value is not None:
                setattr(self, key, roll(value))
        rolled = [[roll(c) for c in Na] for Na in self.N]
        self.N = rolled[-1:] + rolled[:-1]
        self.x = [None if c is None else roll(c) for c in self.x]
        if self.e_up is not None:
            self.e_up = [roll(c) for c in self.e_up]
        self.vaxis, self.haxes = self._shift_axes(+1)

    def _shift_axes(self, step):
        ndim = self.J.ndim
        shift = lambda a: (a + step) % ndim if a is not None else None
        return shift(self.vaxis), tuple(shift(a) for a in self.haxes)


def apply_gradient_map(metric, dp):
    """Physical gradients from computational ones: dp <- A @ dp (in place).

    General curvilinear chain rule: with N_a = J grad xi_a (the duality
    identity, Phase 1) the physical gradient is

        (grad_x p)_k = sum_a A_{k a} dp/dxi_a,    A_{k a} = (N_a)_k / J.

    ``dp`` is the axis-indexed list of the cell-gradient arrays (entries
    beyond ndim — quasi-2D callers pass three — are untouched); on return
    entry k carries the Cartesian component k (canonical orientation only,
    like every caller). For the vertical-line metric this reduces to the
    legacy A = [[1, -G1/J, 0], [0, 1/J, 0], [0, -G2/J, 1]] map to within
    one ulp (the diagonal picks up J * (1/J)); with h == 0 it is the exact
    identity. The diagonal term leads each row's contraction so the
    reduction is deterministic on both backends.
    """
    ooJ = metric.ooJ
    N = metric.N
    ndim = metric.J.ndim
    dp_in = [dp[a] for a in range(ndim)]
    for k in range(ndim):
        acc = (N[k][k] * ooJ) * dp_in[k]
        for a in range(ndim):
            if a != k:
                acc = acc + (N[a][k] * ooJ) * dp_in[a]
        dp[k] = acc
    return dp


def _role_axes(metric, nroles):
    """Role -> Cartesian/array axis map (h1, v[, h2]), from the fixed axes."""
    if nroles == 2:
        return (metric.cart_haxes[0], metric.cart_v)
    return (metric.cart_haxes[0], metric.cart_v, metric.cart_haxes[1])


def _fold_normals(metric, h_role, nroles):
    """M = (1/J) N H^-1 N^T, role-indexed: M_rs = ooJ * N_r . H^-1 N_s.

    ``h_role`` is the role-indexed H^-1 (its role components are the
    Cartesian components in role positions); the contraction runs
    role-major with the (r' == r, s' == s) diagonal term first so the
    h == 0 reduction (N_r = e_r) returns ``h_role`` bit-exactly.
    """
    ooJ = metric.ooJ
    ax = _role_axes(metric, nroles)
    n = [[metric.N[ax[r]][ax[c]] for c in range(nroles)] for r in range(nroles)]
    M = [[None] * nroles for _ in range(nroles)]
    for r in range(nroles):
        for s in range(nroles):
            acc = (n[r][r] * h_role[r][s]) * n[s][s]
            for rp in range(nroles):
                for sp in range(nroles):
                    if rp == r and sp == s:
                        continue
                    acc = acc + (n[r][rp] * h_role[rp][sp]) * n[s][sp]
            M[r][s] = ooJ * acc
    return tuple(tuple(row) for row in M)


def elliptic_tensor(metric, h_role):
    """Fold the metric into the role-indexed H^-1 tensor: M = (1/J) N H^-1 N^T.

    Returns M = J A^T H^-1 A (role order h1, v, h2) in its general
    curvilinear form via the duality N_a = J grad xi_a: the rhs divergence
    measures D_a(N_a . F) and the momentum correction applies
    H^-1 A grad p, so their composition carries exactly this tensor. For
    the vertical-line metric and H^-1 == I it is the classic
    terrain-following tensor

        [[J, -G1, 0], [-G1, (1 + G1^2 + G2^2)/J, -G2], [0, -G2, J]]

    to within one ulp (J (1/J) products), and with h == 0 (N_r = e_r) it
    reduces bit-exactly to ``h_role``.
    """
    return _fold_normals(metric, h_role, 3)


def elliptic_tensor_2d(metric, h2x2):
    """2D restriction of :func:`elliptic_tensor` to the (h1, v) block.

    ``h2x2 = ((h11, h12), (h21, h22))`` is the in-plane H^-1 block (in 2D
    the out-of-plane H^-1 rotation never enters the divergence, so the
    composition carries exactly this 2x2 tensor). Same general fold
    M = (1/J) N H^-1 N^T; with h == 0 it reduces bit-exactly to ``h2x2``.
    """
    return _fold_normals(metric, h2x2, 2)


def elliptic_diag_geometric(metric):
    """Role-ordered diagonal of the geometric tensor M with H^-1 == I.

    M_rr = (1/J) |N_r|^2 — the terrain factors the 2D preconditioner
    folds into its diagonal (H^-1 is deliberately excluded there, so a
    forced-flat metric preconditions bit-identically to the plain path).
    The vertical role reduces to the legacy (1 + G1^2 (+ G2^2)) / J
    bit-exactly; horizontal roles give (1/J) J^2 = J to within one ulp.
    """
    ooJ = metric.ooJ
    ndim = metric.J.ndim
    ax = _role_axes(metric, ndim)
    out = []
    for r in range(ndim):
        comps = metric.N[ax[r]]
        acc = comps[ax[r]] ** 2
        for c in range(ndim):
            if c != ax[r]:
                acc = acc + comps[c] ** 2
        out.append(ooJ * acc)
    return tuple(out)


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


def _effective_callable(ud, a_h1, a_h2, fn):
    """``fn`` with periodic-wrapped arguments (h1, h2 role order)."""
    wrap1 = _coordinate_wrap(ud, a_h1)
    wrap2 = _coordinate_wrap(ud, a_h2)
    return lambda xi1, xi2: fn(wrap1(xi1), wrap2(xi2))


def _effective_orography(ud, a_h1, a_h2):
    """ud.orography with periodic-wrapped arguments (h1, h2 role order)."""
    return _effective_callable(ud, a_h1, a_h2, ud.orography)


def _coord_view(grid_obj, axis, ndim):
    """Coordinate array of `axis`, shaped to broadcast over an ndim field."""
    shape = [1] * ndim
    shape[axis] = -1
    return axes.coords_along(grid_obj, axis).reshape(shape)


def _terrain_slope(ud, heff, grad, a_h1, a_h2, xi1, xi2, which, spacing):
    """dh/dxi_which (role index 0 or 1): analytic if provided, else FD.

    ``grad`` is the role-ordered tuple of analytic gradient callables for
    the orography ``heff`` wraps (or None for FD). Both paths wrap
    periodic coordinates: the analytic gradient is evaluated at the
    wrapped points, the central difference differentiates the wrapped
    (periodic) effective orography so the seam is consistent.
    """
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

    grad = getattr(ud, "orography_grad", None)

    def slope_of(which, spacing):
        return _terrain_slope(ud, heff, grad, a_h1, a_h2, xi1, xi2, which, spacing)

    n_comp = getattr(transform, "n_components", 1)
    if n_comp == 2:
        # two-scale transforms (SLEVE): h splits into a smooth part and the
        # residual. The smooth part MUST be wrapped with the same coordinate
        # map as the total, or the residual breaks the periodic-seam
        # consistency of the metric. No silent default: a missing split
        # would degenerate SLEVE into something else entirely.
        smooth = getattr(ud, "orography_smooth", None)
        if smooth is None:
            raise ValueError(
                "two-component vertical transform requires ud.orography_smooth "
                "(the large-scale part of ud.orography; the residual is the "
                "small-scale component)"
            )
        heff_s = _effective_callable(ud, a_h1, a_h2, smooth)
        h_s = heff_s(xi1, xi2)
        grad_s = getattr(ud, "orography_smooth_grad", None)

        h_arg = (h_s, h - h_s)

        def dh_of(which, spacing):
            dh_tot = slope_of(which, spacing)
            dh_s = _terrain_slope(
                ud, heff_s, grad_s, a_h1, a_h2, xi1, xi2, which, spacing
            )
            return (dh_s, dh_tot - dh_s)

    else:
        h_arg = h
        dh_of = slope_of

    def full(expr):
        return np.ascontiguousarray(
            np.broadcast_to(expr, shape).astype(np.float64, copy=False)
        )

    J = full(transform.jacobian(eta, h_arg, eta0, etat))
    if np.any(J <= 0.0):
        raise ValueError(
            "terrain transform produced non-positive Jacobian: "
            "orography reaches or exceeds the domain top"
        )

    z = full(transform.z(eta, h_arg, eta0, etat))

    G1 = full(transform.slope(eta, dh_of(0, grid_obj.dxyz[a_h1]), eta0, etat))
    if a_h2 is not None:
        G2 = full(transform.slope(eta, dh_of(1, grid_obj.dxyz[a_h2]), eta0, etat))
    else:
        G2 = None

    return MetricFields(J=J, G1=G1, G2=G2, z=z, vaxis=v, haxes=(a_h1, a_h2))


class CurvilinearMap:
    """Analytic curvilinear map x(xi) for the general metric path.

    Tier 1/2 of the Klein generalization: maps that keep gravity a
    coordinate direction, x = (x(xi_0), ..., z(..., eta, ...), ...).
    Subclasses implement elementwise, broadcastable methods of the
    computational coordinates ``xi`` (a list indexed by ARRAY AXIS in
    canonical orientation, i.e. Cartesian order; entries broadcast over
    the grid like the builder's coordinate views).

    ``vertical_line`` (class or instance attribute, default True)
    declares whether the map's vertical coordinate lines are parallel
    Cartesian lines. Tier-3 maps (sphere) set it False: the builder then
    leaves the legacy slope scalars G1/G2 unset (undefined for such maps)
    and derives the up-direction ``e_up`` from the vertical normal.
    """

    vertical_line = True

    def coordinates(self, xi):
        """Physical coordinates x_k(xi) as a Cartesian-indexed list;
        entries may be ``None`` where the map is the identity."""
        raise NotImplementedError

    def tangents(self, xi):
        """Tangents t_a = dx/dxi_a: a list over computational axes ``a``
        of Cartesian-component lists (each entry broadcastable)."""
        raise NotImplementedError

    def height(self, xi):
        """Generalized altitude (gravity potential coordinate) as one
        broadcastable expression, or ``None`` to default to the vertical
        physical coordinate x[v]. A spherical map returns ``r - a``."""
        return None

    def up_direction(self, xi):
        """GRAVITY direction as a Cartesian-component list, or ``None``
        to default to the vertical coordinate-surface normal N_v/|N_v|.
        Terrain maps must override: their gravity stays radial while the
        coordinate surfaces tilt with the slope."""
        return None


def _cross_components(a, b):
    """Elementwise cross product of Cartesian-component triples."""
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def build_metric_fields_from_map(grid_obj, ud, cmap):
    """Build MetricFields for one grid from an analytic CurvilinearMap.

    The general-path sibling of :func:`build_metric_fields`: area normals
    N_a from cross products of the map tangents, J from the triple
    product, physical coordinates from the map. The legacy scalars are
    the EFFECTIVE values the unconverted consumers need —
    z = x[v], G_h = -(N_v)_h / (N_v)_v (exact slopes for vertical-line
    maps, effective slopes under horizontal stretching) — so the wall
    reflection and hydrostate sampling stay correct through Tier 2.

    For a vertical-line map this reproduces :func:`build_metric_fields`
    bit-exactly (the Phase-0 reduction contract).
    """
    ndim = grid_obj.ndim
    v = axes.vertical_axis(ud)
    if ndim == 2:
        a_h1, a_h2 = 0, None
    else:
        a_h1, a_h2 = axes.horizontal_axes(v)

    shape = tuple(int(grid_obj.sc[dim]) for dim in range(ndim))
    xi = [_coord_view(grid_obj, a, ndim) for a in range(ndim)]

    def full(expr):
        return np.ascontiguousarray(
            np.broadcast_to(expr, shape).astype(np.float64, copy=False)
        )

    t = cmap.tangents(xi)
    if ndim == 2:
        # 2D normals: N_a = J grad xi_a is the (rotated) other tangent
        J = full(t[0][0] * t[1][1] - t[0][1] * t[1][0])
        N = [
            [full(t[1][1]), full(-t[1][0])],
            [full(-t[0][1]), full(t[0][0])],
        ]
    else:
        N_raw = [
            _cross_components(t[1], t[2]),
            _cross_components(t[2], t[0]),
            _cross_components(t[0], t[1]),
        ]
        J = full(t[0][0] * N_raw[0][0] + t[0][1] * N_raw[0][1] + t[0][2] * N_raw[0][2])
        N = [[full(c) for c in Na] for Na in N_raw]
    if np.any(J <= 0.0):
        raise ValueError(
            "curvilinear map produced non-positive Jacobian: "
            "the map must be orientation-preserving everywhere"
        )

    x = [None if c is None else full(c) for c in cmap.coordinates(xi)]
    z = x[v] if x[v] is not None else full(xi[v] + 0.0 * J)

    vertical_line = getattr(cmap, "vertical_line", True)
    if vertical_line:
        # effective slopes off the vertical normal (see docstring)
        G1 = full(-N[v][a_h1] / N[v][v])
        G2 = full(-N[v][a_h2] / N[v][v]) if a_h2 is not None else None
        e_up = None
        # t_v = (0, dz/deta, 0): |t_v| is exactly the v-component (no
        # sqrt, keeps the Phase-0 reduction bit-exact)
        h_v = full(t[v][v] + 0.0 * J)
    else:
        # slope scalars are undefined ((N_v)_k passes through zero); the
        # up direction is carried as data instead. The GRAVITY direction
        # is the map's up_direction hook when provided (a terrain map's
        # gravity stays radial while its coordinate-surface normal tilts
        # with the slope); the surface normal N_v/|N_v| is the default
        # (they coincide for terrain-free maps).
        G1 = None
        G2 = None
        up_expr = cmap.up_direction(xi) if hasattr(cmap, "up_direction") else None
        if up_expr is not None:
            e_up = [full(c) for c in up_expr]
        else:
            norm_v = np.sqrt(sum(np.asarray(c) ** 2 for c in N[v]))
            e_up = [full(np.asarray(c) / norm_v) for c in N[v]]
        h_v = full(np.sqrt(sum(np.asarray(c) ** 2 for c in t[v])))

    height_expr = cmap.height(xi) if hasattr(cmap, "height") else None
    height = z if height_expr is None else full(height_expr)

    return MetricFields(
        J=J,
        G1=G1,
        G2=G2,
        z=z,
        vaxis=v,
        haxes=(a_h1, a_h2),
        N=N,
        x=x,
        height=height,
        h_v=h_v,
        e_up=e_up,
        vertical_line=vertical_line,
    )
