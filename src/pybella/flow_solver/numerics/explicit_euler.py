import numpy as np

from ...utils.operators import convolution, divergence, gradient
from ..utils.boundary import cell_boundary as bdry_c
from ..utils.boundary import node_boundary as bdry_n
from ..utils.boundary import common as bdry


def do_forward_step(mem, ud, dt, writer=None, label=None, debug=False):
    # Unpack frequently used variables
    th, sol, npf, node, elem = mem.th, mem.sol, mem.npf, mem.node, mem.elem
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
    p2n = npf.p2_nodes
    dp2n = np.zeros_like(p2n)

    S0c = npf.HydroState.get_S0c(elem)
    dSdy = npf.HydroState_n.get_dSdy(elem, node)

    # Compute divergence
    npf.rhs[...] = divergence.compute_at_nodes(npf.rhs, elem, sol, ud)
    if not hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        bdry.scale_wall_node_values(npf.rhs, node, ud, 2.0)

    if debug:
        writer.populate(str(label), "rhs", npf.rhs)

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

    # the w-row applies in 2D too: dpdz is zero there, but the Coriolis terms
    # are not. Restricting it to ndim == 3 gave 2D runs only the implicit
    # half of the out-of-plane Coriolis rotation — found 2026-06-09 by the
    # Baldauf-Brdar analytic oracle (w_out error pinned at ~0.44 rel-L2 with
    # a sim/ref amplitude ratio ~0.6, independent of dt).
    rhow -= dt * (rhoYovG * dpdz - corr_v * drhou + corr_h1 * drhov)

    # Scalar update (rhoX)
    sol.rhoX[...] = (rho * (rho / rhoY - S0c)) - dt * (v * dSdy) * rho

    # Compressibility correction to p2
    dp2n[node.i1] -= dt * dpidP * npf.rhs
    npf.p2_nodes += ud.compressibility * dp2n

    # Boundary conditions
    bdry_n.set_ghost_nodes(npf.p2_nodes, node, ud)
