import numba as nb

from ...utils import axes


# Refactored main function
def multiply_inverse_terms(
    Vec, mem, ud, dt, attrs=("rhou", "rhov", "rhow"), get_coeffs=False
):
    """Apply H^-1 (Coriolis/buoyancy coupling matrix inverse) to a vector field.

    ``attrs`` is AXIS-indexed (u, v, w component names); the role binding
    (h1, vertical, h2) happens here via the cyclic axis permutation, so the
    njit kernels below — written in role symbols (wh1, wv, wh2) — stay
    unchanged for any vertical axis.
    """
    if getattr(ud, "backend", "numpy") == "jax":
        from ...backends.jax_ops import coriolis as jax_coriolis

        return jax_coriolis.multiply_inverse_terms(
            Vec, mem, ud, dt, attrs=attrs, get_coeffs=get_coeffs
        )

    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[axes.vertical_axis(ud)]
    Msq = ud.Msq

    ax_h1, ax_v, ax_h2 = axes.role_perm(axes.vertical_axis(ud))
    wdt = dt * ud.coriolis_strength
    wh1, wv, wh2 = wdt[ax_h1], wdt[ax_v], wdt[ax_h2]
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
    if getattr(ud, "backend", "numpy") == "jax":
        from ...backends.jax_ops import coriolis as jax_coriolis

        return jax_coriolis.compute_inverse_coefficients(mem, ud, dt)

    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[axes.vertical_axis(ud)]
    Msq = ud.Msq

    ax_h1, ax_v, ax_h2 = axes.role_perm(axes.vertical_axis(ud))
    wdt = dt * ud.coriolis_strength
    wh1, wv, wh2 = wdt[ax_h1], wdt[ax_v], wdt[ax_h2]
    strat = mem.npf.HydroState_n.get_dSdy(mem.elem, mem.node)
    Y = mem.sol.rhoY / mem.sol.rho
    nu = -(dt**2) * (g / Msq) * strat * Y

    views = mem.cache.get_coriolis_array_views(nu.shape)
    h11, h12, h13, h21, h22, h23, h31, h32, h33, denom = views
    _compute_coriolis_coefficients(
        h11, h12, h13, h21, h22, h23, h31, h32, h33, denom, wh1, wh2, wv, nu, nonhydro
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
    h22[...] = nonhydro * (1.0 + wv_sq) * denom
    h23[...] = (wh2 * wv + wh1) * denom

    # Row 3: W equation coefficients
    h31[...] = (wh1 * wh2 + nu_nh * wv) * denom
    h32[...] = nonhydro * (wh2 * wv - wh1) * denom
    h33[...] = (nu_nh + wh2_sq) * denom

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
