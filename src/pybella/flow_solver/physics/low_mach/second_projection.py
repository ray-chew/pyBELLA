import itertools as it

import numpy as np
import scipy as sp
from numba import njit

from ....utils import options as opts
from ....utils import operators

from ...utils import boundary as bdry
from . import laplacian as lm_lp


class solver_counter(object):
    """
    taken from https://stackoverflow.com/questions/33512081/getting-the-number-of-iterations-of-scipys-gmres-iterative-method

    """

    def __init__(self, disp=True):
        self.niter = 0

    def __call__(self, rk=None):
        self.niter += 1
        self.rk = rk


def euler_forward_non_advective(mem, ud, dt, writer=None, label=None, debug=False):
    # Unpack frequently used variables
    th, sol, mpv, node, elem = mem.th, mem.sol, mem.mpv, mem.node, mem.elem
    ndim = elem.ndim

    nonhydro = ud.nonhydrostasy
    g, Msq = ud.gravity_strength[1], ud.Msq
    Ginv = th.Gammainv
    corr_h1, corr_v, corr_h2 = ud.coriolis_strength
    u0, v0, w0 = ud.u_wind_speed, ud.v_wind_speed, ud.w_wind_speed

    # Reusable derived quantities
    rho, rhoY, rhoX = sol.rho, sol.rhoY, sol.rhoX
    rhou, rhov, rhow = sol.rhou, sol.rhov, sol.rhow

    # Pressure and derivatives
    p2n = mpv.p2_nodes
    dp2n = np.zeros_like(p2n)

    S0c = mpv.HydroState.get_S0c(elem)
    dSdy = mpv.HydroState_n.get_dSdy(elem, node)

    # Compute divergence
    mpv.rhs[...] = divergence_nodes(mpv.rhs, elem, sol, ud)
    if not hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        bdry.scale_wall_node_values(mpv.rhs, node, ud, 2.0)

    if debug:
        writer.populate(str(label), "rhs", mpv.rhs)

    # Compute compressibility kernel
    kernel = operators.get_averaging_kernel(ndim, width=2)
    dpidP = (th.gm1 / Msq) * operators.apply_convolution_kernel(
        rhoY ** (th.gamm - 2.0), kernel=kernel, normalize=True, use_numba=True
    )

    rhoYovG = Ginv * rhoY
    dbuoy = rhoY * (rhoX / rho)

    # Pressure gradients
    dpdx, dpdy, dpdz = operators.compute_gradient_nodes(p2n, ndim, node.dxyz)

    # Wind perturbations
    drhou = rhou - u0 * rho
    drhov = rhov - v0 * rho
    drhow = rhow - w0 * rho
    v = rhov / rho

    # Momentum update (u, v, w)
    rhou -= dt * (rhoYovG * dpdx - corr_h2 * drhov + corr_v * drhow)
    rhov -= (
        dt
        * (
            rhoYovG * dpdy
            + (g / Msq) * dbuoy * nonhydro
            - corr_h1 * drhow
            + corr_h2 * drhou
        )
        * (1 - ud.is_ArakawaKonor)
    )

    if ndim == 3:
        rhow -= dt * (rhoYovG * dpdz - corr_v * drhou + corr_h1 * drhov)

    # Scalar update (rhoX)
    sol.rhoX[...] = (rho * (rho / rhoY - S0c)) - dt * (v * dSdy) * rho

    # Compressibility correction to p2
    dp2n[node.i1] -= dt * dpidP * mpv.rhs
    mpv.p2_nodes += ud.compressibility * dp2n

    # Boundary conditions
    bdry.set_ghostnodes_p2(mpv.p2_nodes, node, ud)
    bdry.set_explicit_boundary_data(sol, elem, ud, th, mpv)


def euler_backward_non_advective_expl_part(mem, ud, dt):
    nonhydro = ud.nonhydrostasy
    g = ud.gravity_strength[1]
    Msq = ud.Msq

    dbuoy = mem.sol.rhoY * (mem.sol.rhoX / mem.sol.rho)
    mem.sol.rhov = (nonhydro * mem.sol.rhov) - dt * (g / Msq) * dbuoy

    mem.sol.mod_bg_wind(ud, -1.0)

    multiply_inverse_coriolis(mem.sol, mem, ud, dt)

    mem.sol.mod_bg_wind(ud, +1.0)

    bdry.set_explicit_boundary_data(mem.sol, mem.elem, ud, mem.th, mem.mpv)


def euler_backward_non_advective_impl_part(
    mem,
    ud,
    dt,
    Sol0=None,
    writer=None,
    label=None,
    debug=False,
):
    if not debug:
        writer = None
    nc = mem.node.sc

    if writer != None:
        writer.populate(str(label), "p2_initial", mem.mpv.p2_nodes)

    if Sol0 is not None:
        bdry.set_explicit_boundary_data(Sol0, mem.elem, ud, mem.th, mem.mpv)
        operator_coefficients_nodes(mem, ud, dt)
    else:
        bdry.set_explicit_boundary_data(mem.sol, mem.elem, ud, mem.th, mem.mpv)
        operator_coefficients_nodes(mem, ud, dt)

    if writer != None:
        writer.populate(str(label), "hcenter", mem.mpv.wcenter)
        writer.populate(str(label), "wplusx", mem.mpv.wplus[0])
        writer.populate(str(label), "wplusy", mem.mpv.wplus[1])
        (
            writer.populate(str(label), "wplusz", mem.mpv.wplus[2])
            if mem.elem.ndim == 3
            else writer.populate(str(label), "wplusz", np.zeros_like(mem.mpv.wplus[0]))
        )

    bdry.set_ghostnodes_p2(mem.mpv.p2_nodes, mem.node, ud)
    correction_nodes(mem, ud, dt, mem.mpv.p2_nodes, 0)
    bdry.set_explicit_boundary_data(mem.sol, mem.elem, ud, mem.th, mem.mpv)

    mem.mpv.rhs[...] = divergence_nodes(mem.mpv.rhs, mem.elem, mem.sol, ud)

    if writer != None:
        writer.populate(str(label), "rhs", mem.mpv.rhs)

    mem.mpv.rhs /= dt

    if ud.is_compressible == 0:
        if ud.is_ArakawaKonor:
            mem.mpv.rhs -= mem.mpv.wcenter * mem.mpv.dp2_nodes
            mem.mpv.wcenter[...] = 0.0
        else:
            mem.mpv.rhs = (
                ud.compressibility * mem.mpv.rhs
                + (1.0 - ud.compressibility) * mem.mpv.rhs
            )
            mem.mpv.wcenter[...] *= ud.compressibility
    else:
        mem.mpv.wcenter *= ud.compressibility

    if writer != None:
        writer.populate(str(label), "rhs_nodes", mem.mpv.rhs)

    mem.mpv.rhs[...] = mem.mpv.rhs

    # prepare initial left-hand side and the laplacian stencil
    if mem.elem.ndim == 2:
        Vec = mem.mpv
        coriolis_params = multiply_inverse_coriolis(
            Vec, mem, ud, dt, attrs=("u", "v", "w"), get_coeffs=True
        )

        diag_inv = lm_lp.precon_diag_prepare(mem.mpv, mem.node)
        mem.mpv.rhs *= diag_inv

        p2 = mem.mpv.p2_nodes[mem.node.i2].T
        lap = lm_lp.get_lap2D(mem.mpv, mem.node, coriolis_params, diag_inv, ud)
        sh = p2.shape[0] * p2.shape[1]

    elif mem.elem.ndim == 3:
        lap = lm_lp.stencil_27pt(mem.elem, mem.node, mem.mpv, ud, diag_inv, dt)
        sh = p2.reshape(-1).shape[0]

    lap = sp.sparse.linalg.LinearOperator((sh, sh), lap)
    # lap = LinearOperator(sh,lap)

    counter = solver_counter()

    # prepare right-hand side
    if mem.elem.ndim == 2:
        rhs_inner = mem.mpv.rhs[mem.node.i1].T.ravel()
    else:
        rhs_inner = mem.mpv.rhs[mem.node.i1].ravel()

    p2, _ = sp.sparse.linalg.bicgstab(
        lap, rhs_inner, atol=ud.tol, maxiter=ud.max_iterations, callback=counter
    )

    p2_full = np.zeros(nc).squeeze()
    if mem.elem.ndim == 2:
        p2_full[mem.node.i2] = p2.reshape(mem.mpv.rhs[mem.node.i1].T.shape).T
    elif mem.elem.ndim == 3:
        p2_full[mem.node.i1] = p2.reshape(ud.inx + 2, ud.iny + 2, ud.inz + 2)

    if writer != None:
        writer.populate(str(label), "p2_full", p2_full)

    bdry.set_ghostnodes_p2(p2_full, mem.node, ud)
    correction_nodes(mem, ud, dt, p2_full, 1)

    mem.mpv.p2_nodes[...] += p2_full
    bdry.set_ghostnodes_p2(mem.mpv.p2_nodes, mem.node, ud)
    bdry.set_explicit_boundary_data(mem.sol, mem.elem, ud, mem.th, mem.mpv)


def correction_nodes(mem, ud, dt, p, updt_chi):
    ndim = mem.node.ndim
    Gammainv = mem.th.Gammainv

    dSdy = mem.mpv.HydroState_n.get_dSdy(mem.elem, mem.node)

    Dpx, Dpy, Dpz = operators.compute_gradient_nodes(p, mem.elem.ndim, mem.node.dxyz)

    thinv = mem.sol.rho / mem.sol.rhoY

    Y = mem.sol.rhoY / mem.sol.rho
    coeff = Gammainv * mem.sol.rhoY * Y

    mem.mpv.u[...] = -dt * coeff * Dpx
    mem.mpv.v[...] = -dt * coeff * Dpy
    mem.mpv.w[...] = -dt * coeff * Dpz

    multiply_inverse_coriolis(mem.mpv, mem, ud, dt, attrs=["u", "v", "w"])

    mem.sol.rhou += thinv * mem.mpv.u
    mem.sol.rhov += thinv * mem.mpv.v
    mem.sol.rhow += thinv * mem.mpv.w if ndim == 3 else 0.0
    mem.sol.rhoX += -updt_chi * dt * dSdy * mem.sol.rhov

    assert True


def operator_coefficients_nodes(mem, ud, dt):
    Gammainv = mem.th.Gammainv
    ndim = mem.node.ndim

    ccenter = -ud.Msq * mem.th.gm1inv / (dt**2)
    cexp = 2.0 - mem.th.gamm

    Y = mem.sol.rhoY / mem.sol.rho
    coeff = Gammainv * mem.sol.rhoY * Y

    for dim in range(ndim):
        mem.mpv.wplus[dim][...] = coeff

    kernel = operators.get_averaging_kernel(ndim, width=2)

    mem.mpv.wcenter = ccenter * operators.apply_convolution_kernel(
        mem.sol.rhoY**cexp, kernel
    )

    if not hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        bdry.scale_wall_node_values(mem.mpv.wcenter, mem.node, ud)


@njit(cache=True)
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


@njit(cache=True)
def apply_coriolis_matrix_inplace(
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


# Refactored main function
def multiply_inverse_coriolis(
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

    apply_coriolis_matrix_inplace(VecU, VecV, VecW, U, V, W, wh1, wh2, wv, nu, nonhydro)

    # Return coefficients
    if get_coeffs:
        h11, h12, _, h21, h22, h23, h31, h32, h33 = _compute_coriolis_coefficients(
            wh1, wh2, wv, nu, nonhydro
        )
        return (h11.T, h22.T, h12.T, h21.T)


def divergence_nodes(rhs, elem, sol, ud):
    """Main divergence function - handles boundary conditions and calls JIT-compiled core."""
    ndim = elem.ndim

    # Handle boundary conditions
    if not hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        if (
            ud.bdry_type[1] == opts.BdryType.WALL
            or ud.bdry_type[1] == opts.BdryType.RAYLEIGH
        ):
            sol.rhou[:, :2, ...] = 0.0
            sol.rhov[:, :2, ...] = 0.0
            sol.rhow[:, :2, ...] = 0.0
            sol.rhou[:, -2:, ...] = 0.0
            sol.rhov[:, -2:, ...] = 0.0
            sol.rhow[:, -2:, ...] = 0.0

    # Call appropriate JIT-compiled function
    if ndim == 2:
        rhs[:] = _momentum_pot_temp_divergence_2d_jit(
            sol.rho, sol.rhou, sol.rhov, sol.rhoY, elem.dx, elem.dy
        )
    else:
        _momentum_pot_temp_divergence_3d_jit(
            rhs,
            sol.rho,
            sol.rhou,
            sol.rhov,
            sol.rhow,
            sol.rhoY,
            elem.dx,
            elem.dy,
            elem.dz,
        )

    return rhs


@njit(cache=True)
def _momentum_pot_temp_divergence_2d_jit(rho, rhou, rhov, rhoY, dx, dy):
    """
    JIT-compiled 2D momentum-potential temperature divergence calculation.
    Computes ∇·(ρu θ, ρv θ) where θ = ρY/ρ is the potential temperature.
    """
    # Calculate potential temperature θ = ρY / ρ
    theta = rhoY / rho

    # Compute momentum-potential temperature flux components
    rhou_theta = rhou * theta  # x-momentum flux weighted by potential temperature
    rhov_theta = rhov * theta  # y-momentum flux weighted by potential temperature

    # Use generic divergence operator
    return operators.compute_divergence_2d(rhou_theta, rhov_theta, dx, dy)


@njit(cache=True)
def _momentum_pot_temp_divergence_3d_jit(rhs, rho, rhou, rhov, rhow, rhoY, dx, dy, dz):
    """
    JIT-compiled 3D momentum-potential temperature divergence calculation.
    Computes ∇·(ρu θ, ρv θ, ρw θ) where θ = ρY/ρ is the potential temperature.
    """
    # Calculate potential temperature θ = ρY / ρ
    theta = rhoY / rho

    # Compute momentum-potential temperature flux components
    rhou_theta = rhou * theta  # x-momentum flux weighted by potential temperature
    rhov_theta = rhov * theta  # y-momentum flux weighted by potential temperature
    rhow_theta = rhow * theta  # z-momentum flux weighted by potential temperature

    # Use generic total divergence operator
    total_div = operators.compute_divergence_3d_total(
        rhou_theta, rhov_theta, rhow_theta, dx, dy, dz
    )

    # Assign to inner region
    rhs[1:-1, 1:-1, 1:-1] = total_div
