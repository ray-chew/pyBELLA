import numba as nb

# Refactored main function
def multiply_inverse_terms(
    Vec, mem, ud, dt, attrs=("rhou", "rhov", "rhow"), get_coeffs=False
):
    """Coriolis matrix multiplication."""
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[1]
    Msq = ud.Msq

    wh1, wv, wh2 = dt * ud.coriolis_strength
    strat = mem.mpv.HydroState_n.get_dSdy(mem.elem, mem.node)
    Y = mem.sol.rhoY / mem.sol.rho
    nu = -(dt**2) * (g / Msq) * strat * Y

    # Get vector components
    VecU = getattr(Vec, attrs[0])
    VecV = getattr(Vec, attrs[1])
    VecW = getattr(Vec, attrs[2])

    U, V, W = mem.cache.get_velocity_array_views(VecU.shape)

    _apply_coriolis_matrix_inplace(VecU, VecV, VecW, U, V, W, wh1, wh2, wv, nu, nonhydro)

    # Return coefficients
    if get_coeffs:
        h11, h12, _, h21, h22, _, _, _, _ = _compute_coriolis_coefficients(
            wh1, wh2, wv, nu, nonhydro
        )
        return (h11.T, h22.T, h12.T, h21.T)
    

@nb.njit(cache=True)
def _compute_coriolis_coefficients(wh1, wh2, wv, nu, nonhydro):
    """Compute coefficients for the H^-1 matrix multiplication.

    This corresponds to equation (C11) in the mathematical formulation.
    """
    # Common terms
    wh1_sq = wh1 * wh1
    wh2_sq = wh2 * wh2
    wv_sq = wv * wv
    nu_nh = nu + nonhydro

    # Denominator (det(H))
    denom = 1.0 / (wh1_sq + wh2_sq + nu_nh * (wv_sq + 1.0))

    # H^-1 matrix elements (row-major order)
    # Row 1: U equation coefficients
    h11 = (wh1_sq + nu_nh) * denom
    h12 = nonhydro * (wh1 * wv + wh2) * denom
    h13 = (wh1 * wh2 - nu_nh * wv) * denom

    # Row 2: V equation coefficients
    h21 = (wh1 * wv - wh2) * denom
    h22 = nonhydro * (1.0 + wv_sq) * denom
    h23 = (wh2 * wv + wh1) * denom

    # Row 3: W equation coefficients
    h31 = (wh1 * wh2 + nu_nh * wv) * denom
    h32 = nonhydro * (wh2 * wv - wh1) * denom
    h33 = (nu_nh + wh2_sq) * denom

    return h11, h12, h13, h21, h22, h23, h31, h32, h33

@nb.njit(cache=True)
def _apply_coriolis_matrix_inplace(
    u_vec, v_vec, w_vec, U, V, W, wh1, wh2, wv, nu, nonhydro
):
    """Apply H^-1 matrix multiplication in-place.

    Corresponds to the equation: U^{n+1} = H^{-1}(U^{n*} - Δt_{cp}(Pθ)^* ∇π^{n+1})
    """
    # Get matrix coefficients
    h11, h12, h13, h21, h22, h23, h31, h32, h33 = _compute_coriolis_coefficients(
        wh1, wh2, wv, nu, nonhydro
    )

    U[...] = u_vec
    V[...] = v_vec
    W[...] = w_vec

    # Matrix multiplication: [U_new, V_new, W_new] = H^-1 @ [U_old, V_old, W_old]
    u_vec[...] = h11 * U + h12 * V + h13 * W
    v_vec[...] = h21 * U + h22 * V + h23 * W
    w_vec[...] = h31 * U + h32 * V + h33 * W