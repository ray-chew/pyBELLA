import numba as nb
import numpy as np

from ...utils import axes
from ...backends import is_jax_backend


def role_components(mem, ud):
    """Role-ordered rotation components (w_h1, w_v, w_h2), UNSCALED by dt.

    Uniform-rotation path (``ud.coriolis_field`` unset): exactly
    ``ud.coriolis_strength`` indexed by the role permutation — three
    scalars.

    Field path (``ud.coriolis_field`` set): a callable of the three
    CARTESIAN coordinates returning the three Cartesian rotation
    components.
    """
    ax_h1, ax_v, ax_h2 = axes.role_perm(axes.vertical_axis(ud))
    field = getattr(ud, "coriolis_field", None)
    if field is None:
        w = ud.coriolis_strength
        return w[ax_h1], w[ax_v], w[ax_h2]

    elem = mem.elem
    cached = getattr(elem, "_coriolis_field_cache", None)
    if cached is None or cached[0] is not field:
        shape = mem.sol.rho.shape
        metric = elem.metric
        if metric is not None and all(c is not None for c in metric.x):
            coords = metric.x
        else:
            # identity map: physical coordinates are the grid coordinates
            coords = []
            for a in range(elem.ndim):
                view_shape = [1] * elem.ndim
                view_shape[a] = -1
                coords.append(axes.coords_along(elem, a).reshape(view_shape))
        w = field(*coords)
        w = [
            np.ascontiguousarray(
                np.broadcast_to(c, shape).astype(np.float64, copy=False)
            )
            for c in w
        ]
        cached = (field, w)
        elem._coriolis_field_cache = cached
    w = cached[1]
    return w[ax_h1], w[ax_v], w[ax_h2]


def _up_role_components(mem, ud):
    """Role-ordered (e_h1, e_v, e_h2) components of the local up direction,
    or None on vertical-line/no-metric runs (up = the role-v axis;
    ``_compute_coriolis_coefficients`` (eq. C11) is exactly that case)."""
    metric = mem.elem.metric
    if metric is None or metric.e_up is None:
        return None
    ax_h1, ax_v, ax_h2 = axes.role_perm(axes.vertical_axis(ud))
    e = metric.e_up
    return e[ax_h1], e[ax_v], e[ax_h2]


def multiply_inverse_terms(
    Vec, mem, ud, dt, attrs=("rhou", "rhov", "rhow"), get_coeffs=False
):
    """Apply H^-1 (Coriolis/buoyancy coupling matrix inverse) to a vector field.

    ``attrs`` is AXIS-indexed (u, v, w component names); the role binding
    (h1, vertical, h2) happens here via the cyclic axis permutation, so the
    njit kernels below — written in role symbols (wh1, wv, wh2) — stay
    unchanged for any vertical axis.
    """
    if is_jax_backend(ud):
        from ...backends.jax_ops import coriolis as jax_coriolis

        return jax_coriolis.multiply_inverse_terms(
            Vec, mem, ud, dt, attrs=attrs, get_coeffs=get_coeffs
        )

    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[axes.vertical_axis(ud)]
    Msq = ud.Msq

    ax_h1, ax_v, ax_h2 = axes.role_perm(axes.vertical_axis(ud))
    w_h1, w_v, w_h2 = role_components(mem, ud)
    wh1, wv, wh2 = dt * w_h1, dt * w_v, dt * w_h2
    strat = mem.npf.HydroState_n.get_dSdy(mem.elem, mem.node)
    Y = mem.sol.rhoY / mem.sol.rho
    nu = -(dt**2) * (g / Msq) * strat * Y
    shp = nu.shape

    # initialise Coriolis cache is not already done
    h11, h12, h13, h21, h22, h23, h31, h32, h33, denom = (
        mem.cache.get_coriolis_array_views(shp)
    )

    # Get vector components by view, in role order
    VecU = getattr(Vec, attrs[ax_h1])
    VecV = getattr(Vec, attrs[ax_v])
    VecW = getattr(Vec, attrs[ax_h2])

    U, V, W = mem.cache.get_velocity_array_views(VecU.shape)

    e_role = _up_role_components(mem, ud)
    if e_role is not None:
        _compute_coriolis_coefficients_general(
            h11,
            h12,
            h13,
            h21,
            h22,
            h23,
            h31,
            h32,
            h33,
            denom,
            wh1,
            wh2,
            wv,
            e_role[0],
            e_role[1],
            e_role[2],
            nu,
            nonhydro,
        )
        U[...] = VecU
        V[...] = VecV
        W[...] = VecW
        VecU[...] = h11 * U + h12 * V + h13 * W
        VecV[...] = h21 * U + h22 * V + h23 * W
        VecW[...] = h31 * U + h32 * V + h33 * W
    else:
        _apply_coriolis_matrix_inplace(
            h11,
            h12,
            h13,
            h21,
            h22,
            h23,
            h31,
            h32,
            h33,
            denom,
            VecU,
            VecV,
            VecW,
            U,
            V,
            W,
            wh1,
            wh2,
            wv,
            nu,
            nonhydro,
        )

    # Return coefficients
    if get_coeffs:
        # 2D-only path (the (h1, v) block); 2D runs force vertical = axis 1
        h11, h12, _, h21, h22, _, _, _, _, _ = mem.cache.get_coriolis_array_views(shp)
        return (h11.T, h22.T, h12.T, h21.T)


def compute_inverse_coefficients(mem, ud, dt):
    """Fill and return the cached role-indexed H^-1 coefficient fields.

    Exactly the coefficients the apply path uses (eq. C11), exposed for the
    full-tensor elliptic operator, which consumes them as stencil
    coefficient fields. Returns the 10 cached views
    (h11, h12, h13, h21, h22, h23, h31, h32, h33, denom), role-indexed,
    shaped like the buoyancy field nu.
    """
    if is_jax_backend(ud):
        from ...backends.jax_ops import coriolis as jax_coriolis

        return jax_coriolis.compute_inverse_coefficients(mem, ud, dt)

    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[axes.vertical_axis(ud)]
    Msq = ud.Msq

    w_h1, w_v, w_h2 = role_components(mem, ud)
    wh1, wv, wh2 = dt * w_h1, dt * w_v, dt * w_h2
    strat = mem.npf.HydroState_n.get_dSdy(mem.elem, mem.node)
    Y = mem.sol.rhoY / mem.sol.rho
    nu = -(dt**2) * (g / Msq) * strat * Y

    views = mem.cache.get_coriolis_array_views(nu.shape)
    h11, h12, h13, h21, h22, h23, h31, h32, h33, denom = views
    e_role = _up_role_components(mem, ud)
    if e_role is not None:
        _compute_coriolis_coefficients_general(
            h11,
            h12,
            h13,
            h21,
            h22,
            h23,
            h31,
            h32,
            h33,
            denom,
            wh1,
            wh2,
            wv,
            e_role[0],
            e_role[1],
            e_role[2],
            nu,
            nonhydro,
        )
    else:
        _compute_coriolis_coefficients(
            h11,
            h12,
            h13,
            h21,
            h22,
            h23,
            h31,
            h32,
            h33,
            denom,
            wh1,
            wh2,
            wv,
            nu,
            nonhydro,
        )
    return views


@nb.njit(cache=True)
def _compute_coriolis_coefficients(
    h11, h12, h13, h21, h22, h23, h31, h32, h33, denom, wh1, wh2, wv, nu, nonhydro
):
    """Compute coefficients for the H^-1 matrix multiplication.

    This corresponds to equation (C11) in the mathematical formulation.
    """
    # Common terms
    wh1_sq = wh1 * wh1
    wh2_sq = wh2 * wh2
    wv_sq = wv * wv
    nu_nh = nu + nonhydro

    # Denominator (det(H))
    denom[...] = 1.0 / (wh1_sq + wh2_sq + nu_nh * (wv_sq + 1.0))

    # H^-1 matrix elements (row-major order)
    # Row 1: U equation coefficients
    h11[...] = (wh1_sq + nu_nh) * denom
    h12[...] = nonhydro * (wh1 * wv + wh2) * denom
    h13[...] = (wh1 * wh2 - nu_nh * wv) * denom

    # Row 2: V equation coefficients
    h21[...] = (wh1 * wv - wh2) * denom
    # h22 is the vertical-momentum self-coupling. It must NOT carry the
    # ``nonhydro`` (alpha_w) factor: alpha_w lives in nu_nh (the denominator
    # inertia) and in the explicit ``nonhydro * vmom`` discard
    # (implicit_euler.do_explicit_part). With the factor, the hydrostatic case
    # (alpha_w = 0) zeroed h22 -> the elliptic operator's vertical coupling
    # cij[v][v] = wplus[v] * h22 vanished, so the hydrostatic pressure solve
    # could not reconstruct a balanced Exner pressure (it only preserved one).
    # Dropping it gives h22 -> (1 + wv^2)/det = 1/nu_nh (no Coriolis), the
    # vertical Laplacian coefficient the thesis hydrostatic balance prescribes
    # (eq. 4.40/4.42). Bit-identical for alpha_w = 1 (the factor was 1 there).
    h22[...] = (1.0 + wv_sq) * denom
    h23[...] = (wh2 * wv + wh1) * denom

    # Row 3: W equation coefficients
    h31[...] = (wh1 * wh2 + nu_nh * wv) * denom
    h32[...] = nonhydro * (wh2 * wv - wh1) * denom
    h33[...] = (nu_nh + wh2_sq) * denom

    return h11, h12, h13, h21, h22, h23, h31, h32, h33


@nb.njit(cache=True)
def _compute_coriolis_coefficients_general(
    h11,
    h12,
    h13,
    h21,
    h22,
    h23,
    h31,
    h32,
    h33,
    denom,
    wh1,
    wh2,
    wv,
    e1,
    e2,
    e3,
    nu,
    nonhydro,
):
    """General (C11) matrix for an arbitrary up-direction e (role comps).

    Two ingredients, both reducing analytically to
    ``_compute_coriolis_coefficients`` (the axis-aligned kernel) for
    e = role-v (gated to <= 1e-14 in test_scripts/test_sphere_metric.py):

    1. H^-1 with H = I + (alpha_w - 1 + nu) e e^T + W_x — the
       buoyancy/hydrostasy rank-one term along the LOCAL up e — via
       Sherman-Morrison on the pure-rotation inverse:

           C^-1 = (I + w w^T - W_x) / (1 + |w|^2),
           H^-1 = C^-1 - mu (C^-1 e)(e^T C^-1) / (1 + mu e^T C^-1 e),

       mu = nu + alpha_w - 1. The axis-aligned kernel's h22 structure (no
       alpha_w factor on the vertical self-coupling) falls out:
       denom_SM - mu * c22 = 1.
    2. The thesis (C11) alpha_w prefactor on the (horizontal <- vertical)
       couplings — in the axis-aligned kernel h12 = alpha_w * (H^-1)_12
       and h32 = alpha_w * (H^-1)_32, so the DIAGNOSTIC up-momentum input
       never feeds the transverse rows in the hydrostatic limit.
       Coordinate-free:

           G = H^-1 - (1 - alpha_w) (I - e e^T) H^-1 (e e^T).

    Role order: w = (wh1, wv, wh2), e = (e1, e2, e3) = (e_h1, e_v, e_h2).
    """
    w1 = wh1
    w2 = wv
    w3 = wh2
    ooD = 1.0 / (1.0 + (w1 * w1 + w2 * w2 + w3 * w3))

    c11 = (1.0 + w1 * w1) * ooD
    c12 = (w1 * w2 + w3) * ooD
    c13 = (w1 * w3 - w2) * ooD
    c21 = (w2 * w1 - w3) * ooD
    c22 = (1.0 + w2 * w2) * ooD
    c23 = (w2 * w3 + w1) * ooD
    c31 = (w3 * w1 + w2) * ooD
    c32 = (w3 * w2 - w1) * ooD
    c33 = (1.0 + w3 * w3) * ooD

    # t = C^-1 e (column), s = e^T C^-1 (row)
    t1 = c11 * e1 + c12 * e2 + c13 * e3
    t2 = c21 * e1 + c22 * e2 + c23 * e3
    t3 = c31 * e1 + c32 * e2 + c33 * e3
    s1 = e1 * c11 + e2 * c21 + e3 * c31
    s2 = e1 * c12 + e2 * c22 + e3 * c32
    s3 = e1 * c13 + e2 * c23 + e3 * c33

    mu = nu + (nonhydro - 1.0)
    denom[...] = 1.0 / (1.0 + mu * (e1 * t1 + e2 * t2 + e3 * t3))

    h11[...] = c11 - mu * t1 * s1 * denom
    h12[...] = c12 - mu * t1 * s2 * denom
    h13[...] = c13 - mu * t1 * s3 * denom
    h21[...] = c21 - mu * t2 * s1 * denom
    h22[...] = c22 - mu * t2 * s2 * denom
    h23[...] = c23 - mu * t2 * s3 * denom
    h31[...] = c31 - mu * t3 * s1 * denom
    h32[...] = c32 - mu * t3 * s2 * denom
    h33[...] = c33 - mu * t3 * s3 * denom

    # (C11) alpha_w prefactor: G = H^-1 - (1-a)(I - e e^T) H^-1 (e e^T);
    # col = H^-1 e, its e-perpendicular part scales the e-input column
    col1 = h11 * e1 + h12 * e2 + h13 * e3
    col2 = h21 * e1 + h22 * e2 + h23 * e3
    col3 = h31 * e1 + h32 * e2 + h33 * e3
    col_par = e1 * col1 + e2 * col2 + e3 * col3
    onema = 1.0 - nonhydro
    cp1 = onema * (col1 - col_par * e1)
    cp2 = onema * (col2 - col_par * e2)
    cp3 = onema * (col3 - col_par * e3)
    h11[...] = h11 - cp1 * e1
    h12[...] = h12 - cp1 * e2
    h13[...] = h13 - cp1 * e3
    h21[...] = h21 - cp2 * e1
    h22[...] = h22 - cp2 * e2
    h23[...] = h23 - cp2 * e3
    h31[...] = h31 - cp3 * e1
    h32[...] = h32 - cp3 * e2
    h33[...] = h33 - cp3 * e3

    return h11, h12, h13, h21, h22, h23, h31, h32, h33


@nb.njit(cache=True)
def _apply_coriolis_matrix_inplace(
    h11,
    h12,
    h13,
    h21,
    h22,
    h23,
    h31,
    h32,
    h33,
    denom,
    u_vec,
    v_vec,
    w_vec,
    U,
    V,
    W,
    wh1,
    wh2,
    wv,
    nu,
    nonhydro,
):
    """Apply H^-1 matrix multiplication in-place.

    Corresponds to the equation: U^{n+1} = H^{-1}(U^{n*} - Δt_{cp}(Pθ)^* ∇π^{n+1})
    """
    # Get matrix coefficients
    _compute_coriolis_coefficients(
        h11, h12, h13, h21, h22, h23, h31, h32, h33, denom, wh1, wh2, wv, nu, nonhydro
    )

    U[...] = u_vec
    V[...] = v_vec
    W[...] = w_vec

    # Matrix multiplication: [U_new, V_new, W_new] = H^-1 @ [U_old, V_old, W_old]
    u_vec[...] = h11 * U + h12 * V + h13 * W
    v_vec[...] = h21 * U + h22 * V + h23 * W
    w_vec[...] = h31 * U + h32 * V + h33 * W
