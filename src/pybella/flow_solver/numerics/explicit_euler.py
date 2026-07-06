import numpy as np

from ...utils import axes
from ...utils.operators import convolution, divergence, gradient
from ..discretisation import terrain
from ..utils.boundary import cell_boundary as bdry_c
from ..utils.boundary import node_boundary as bdry_n
from ..utils.boundary import common as bdry


def do_forward_step(mem, ud, dt, writer=None, label=None, debug=False):
    # Unpack frequently used variables
    th, sol, npf, node, elem = mem.th, mem.sol, mem.npf, mem.node, mem.elem
    ndim = elem.ndim

    nonhydro = ud.nonhydrostasy
    g, Msq = ud.gravity_strength[axes.vertical_axis(ud)], ud.Msq
    Ginv = th.Gammainv
    # role-ordered Coriolis components (h1, v, h2); identity for vertical = 1.
    # Scalars on the legacy path; per-cell fields when ud.coriolis_field is
    # set (spatially varying rotation, e.g. f(phi) e_r on the sphere)
    from . import coriolis as coriolis_mod

    ax_h1, ax_v, ax_h2 = axes.role_perm(axes.vertical_axis(ud))
    corr_h1, corr_v, corr_h2 = coriolis_mod.role_components(mem, ud)
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

    # Pressure gradients (physical: terrain slope/Jacobian correction via A)
    dpdx, dpdy, dpdz = gradient.compute_at_nodes(p2n, ndim, node.dxyz)
    if elem.metric is not None:
        dpdx, dpdy, dpdz = terrain.apply_gradient_map(elem.metric, [dpdx, dpdy, dpdz])

    # Wind perturbations
    drhou = rhou - u0 * rho
    drhov = rhov - v0 * rho
    drhow = rhow - w0 * rho

    # role-ordered views (h1, v, h2) of the axis-indexed component tuples;
    # in-place updates below mutate the underlying sol arrays
    mom = (rhou, rhov, rhow)
    dmom = (drhou, drhov, drhow)
    dpd = (dpdx, dpdy, dpdz)
    mom_h1, mom_v, mom_h2 = mom[ax_h1], mom[ax_v], mom[ax_h2]
    dm_h1, dm_v, dm_h2 = dmom[ax_h1], dmom[ax_v], dmom[ax_h2]
    dp_h1, dp_v, dp_h2 = dpd[ax_h1], dpd[ax_v], dpd[ax_h2]

    e_up = elem.metric.e_up if elem.metric is not None else None

    if e_up is not None:
        # general map (sphere): buoyancy acts along the LOCAL up e_up and
        # the H1b nonhydro factor applies to the e-PARALLEL part of the
        # whole tendency (the faithful generalization of "the vertical
        # row"); the e-perpendicular part is never alpha_w-suppressed.
        # Role-ordered e components (same axis permutation as momenta).
        e_h1, e_v, e_h2 = e_up[ax_h1], e_up[ax_v], e_up[ax_h2]
        # rhoX couples to the PRE-update e-parallel (radial) velocity,
        # mirroring the legacy vel_v = mom_v / rho placement
        vel_up = (mom_h1 * e_h1 + mom_v * e_v + mom_h2 * e_h2) / rho
        buoy = (g / Msq) * dbuoy
        T_h1 = rhoYovG * dp_h1 + buoy * e_h1 - corr_h2 * dm_v + corr_v * dm_h2
        T_v = rhoYovG * dp_v + buoy * e_v - corr_h1 * dm_h2 + corr_h2 * dm_h1
        T_h2 = rhoYovG * dp_h2 + buoy * e_h2 - corr_v * dm_h1 + corr_h1 * dm_v
        T_dot_e = T_h1 * e_h1 + T_v * e_v + T_h2 * e_h2
        fac_par = nonhydro * (1 - ud.is_ArakawaKonor) - 1.0
        mom_h1 -= dt * (T_h1 + fac_par * T_dot_e * e_h1)
        mom_v -= dt * (T_v + fac_par * T_dot_e * e_v)
        mom_h2 -= dt * (T_h2 + fac_par * T_dot_e * e_h2)

        sol.rhoX[...] = (rho * (rho / rhoY - S0c)) - dt * (vel_up * dSdy) * rho

        # Compressibility correction below is branch-shared; skip the
        # legacy rows
        _legacy_rows = False
    else:
        _legacy_rows = True

    if _legacy_rows:
        vel_v = mom_v / rho

        # Momentum update in role space: gravity/buoyancy acts on the vertical
        # row, Coriolis couples the rows pairwise (cross-product structure)
        mom_h1 -= dt * (rhoYovG * dp_h1 - corr_h2 * dm_v + corr_v * dm_h2)
        # The WHOLE vertical-momentum forward update carries the nonhydro
        # (alpha_w) factor — not just the buoyancy. In the hydrostatic limit
        # (alpha_w = 0) the vertical momentum has no prognostic time update at
        # all: it is diagnosed by the implicit hydrostatic balance solve.
        # Applying the vertical pressure-gradient kick here for alpha_w = 0
        # (as the pre-2026 refactor did) injects a spurious kick on the
        # diagnosed w every corrector, seeding a 2*dt computational mode.
        # Matches the thesis-era euler_forward_non_advective (rhov update *
        # nonhydro). Bit-identical for alpha_w = 1.
        # See dev_notes/hydrostatic_blending.md (Phase H1b).
        mom_v -= (
            dt
            * (rhoYovG * dp_v + (g / Msq) * dbuoy - corr_h1 * dm_h2 + corr_h2 * dm_h1)
            * nonhydro
            * (1 - ud.is_ArakawaKonor)
        )

        # the h2-row applies in 2D too: its pressure gradient is zero there,
        # but the Coriolis terms are not. Restricting it to ndim == 3 gave 2D
        # runs only the implicit half of the out-of-plane Coriolis rotation —
        # found 2026-06-09 by the Baldauf-Brdar analytic oracle (w_out error
        # pinned at ~0.44 rel-L2 with a sim/ref amplitude ratio ~0.6,
        # independent of dt).
        mom_h2 -= dt * (rhoYovG * dp_h2 - corr_v * dm_h1 + corr_h1 * dm_v)

        # Scalar update (rhoX): stratification couples to the vertical velocity
        sol.rhoX[...] = (rho * (rho / rhoY - S0c)) - dt * (vel_v * dSdy) * rho

    # Compressibility correction to p2; with terrain npf.rhs carries J*div F,
    # so the pointwise pi update needs the plain divergence back (1/J_n)
    if node.metric is not None:
        dp2n[node.i1] -= dt * dpidP * (npf.rhs * node.metric.ooJ[node.i1])
    else:
        dp2n[node.i1] -= dt * dpidP * npf.rhs
    npf.p2_nodes += ud.compressibility * dp2n

    # Boundary conditions
    bdry_n.set_ghost_nodes(npf.p2_nodes, node, ud)
