import numpy as np

from ...utils.operators import convolution, divergence, gradient
from ..utils import boundary as bdry

def do_forward_step(mem, ud, dt, writer=None, label=None, debug=False):
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
    mpv.rhs[...] = divergence.compute_at_nodes(mpv.rhs, elem, sol, ud)
    if not hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        bdry.scale_wall_node_values(mpv.rhs, node, ud, 2.0)

    if debug:
        writer.populate(str(label), "rhs", mpv.rhs)

    # Compute compressibility kernel
    kernel = convolution.get_averaging_kernel(ndim, width=2)
    dpidP = (th.gm1 / Msq) * convolution.apply_convolution_kernel(
        rhoY ** (th.gamm - 2.0), kernel=kernel, normalize=True, use_numba=True
    )

    rhoYovG = Ginv * rhoY
    dbuoy = rhoY * (rhoX / rho)

    # Pressure gradients
    dpdx, dpdy, dpdz = gradient.compute_at_nodes(p2n, ndim, node.dxyz)

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
