import numpy as np
import numba as nb

from ...utils import axes
from ..utils.boundary import node_boundary as bdry_n

# Fine auxiliary 1D z-grid for the terrain hydrostate quadrature: at least this
# many points, and at least this many per vertical cell, so the trapezoidal
# integral of 1/stratification converges well below the solver tolerance.
_HYDRO_QUAD_MIN_POINTS = 2048
_HYDRO_QUAD_POINTS_PER_CELL = 16


def column(HydroState, HydroState_n, Y, Y_n, elem, node, th, ud):
    """2D x-y initial-condition helper (vertical = axis 1 by convention)."""
    assert elem.ndim == 2, "column() is a 2D x-y IC helper"
    assert elem.metric is None, "column() does not support terrain"
    Gamma = th.gm1 / th.gamm
    gamm = th.gamm
    gm1 = th.gm1
    Gamma_inv = 1.0 / Gamma
    gm1_inv = 1.0 / gm1

    icy = elem.icy
    igy = elem.igy

    xc_idx = slice(0, -1)
    yc_idx = slice(0, -1)

    c_idx = (xc_idx, yc_idx)

    rhoY0 = 1.0

    g = ud.gravity_strength[axes.vertical_axis(ud)]

    p0 = rhoY0**gamm
    pi0 = rhoY0**gm1
    HydroState_n.rho0[xc_idx, igy] = rhoY0 / Y_n[:, igy]
    HydroState_n.rhoY0[xc_idx, igy] = rhoY0
    HydroState_n.Y0[xc_idx, igy] = Y_n[:, igy]
    HydroState_n.S0[xc_idx, igy] = 1.0 / Y_n[:, igy]
    HydroState_n.p0[xc_idx, igy] = p0
    HydroState_n.p20[xc_idx, igy] = pi0 / ud.Msq

    dys = np.array(
        [-elem.dy] + [-elem.dy / 2] + [elem.dy / 2] + list(np.ones((icy - 3)) * elem.dy)
    )
    S_p = 1.0 / Y[:, :]
    S_m = np.zeros_like(S_p)
    S_m[:, igy - 1 : igy + 1] = 1.0 / Y_n[:, igy].reshape(-1, 1)
    S_m[:, 0] = 1.0 / Y[:, igy - 1]
    S_m[:, igy + 1 :] = 1.0 / Y[:, igy:-1]

    S_integral_p = dys * 0.5 * (S_p + S_m)
    S_integral_p[:, :igy] = np.cumsum(S_integral_p[:, :igy][:, ::-1], axis=1)[:, ::-1]
    S_integral_p[:, igy:] = np.cumsum(S_integral_p[:, igy:], axis=1)

    pi_hydro = pi0 - Gamma * g * S_integral_p
    p_hydro = pi_hydro**Gamma_inv
    rhoY_hydro = pi_hydro**gm1_inv

    HydroState.rho0[c_idx] = rhoY_hydro * S_p
    HydroState.p0[c_idx] = p_hydro
    HydroState.p20[c_idx] = pi_hydro / ud.Msq
    HydroState.S0[c_idx] = S_p
    HydroState.S10[c_idx] = 0.0
    HydroState.Y0[c_idx] = 1.0 / S_p
    HydroState.rhoY0[c_idx] = rhoY_hydro

    Sn_p = 1.0 / Y[:, :]
    dys = np.ones((icy)) * elem.dy
    dys[:igy] *= -1
    Sn_integral_p = dys * Sn_p
    Sn_integral_p[:, :igy] = np.cumsum(Sn_integral_p[:, :igy][:, ::-1], axis=1)[:, ::-1]
    Sn_integral_p[:, igy:] = np.cumsum(Sn_integral_p[:, igy:], axis=1)

    pi_hydro_n = pi0 - Gamma * g * Sn_integral_p
    rhoY_hydro_n = pi_hydro_n**gm1_inv

    HydroState_n.rhoY0[xc_idx, :igy] = rhoY_hydro_n[:, :igy]
    HydroState_n.Y0[xc_idx, :igy] = Y_n[0, :igy]
    HydroState_n.S0[xc_idx, :igy] = 1.0 / Y_n[:, :igy]
    HydroState_n.p0[xc_idx, :igy] = rhoY_hydro_n[:, :igy] ** th.gamm
    HydroState_n.p20[xc_idx, :igy] = pi_hydro_n[:, :igy] / ud.Msq

    HydroState_n.rhoY0[xc_idx, igy + 1 :] = rhoY_hydro_n[:, igy:]
    HydroState_n.Y0[xc_idx, igy + 1 :] = Y_n[0, igy:]
    HydroState_n.S0[xc_idx, igy + 1 :] = 1.0 / Y_n[:, igy:]
    HydroState_n.p0[xc_idx, igy + 1 :] = rhoY_hydro_n[:, igy:] ** th.gamm
    HydroState_n.p20[xc_idx, igy + 1 :] = pi_hydro_n[:, igy:] / ud.Msq


def integrated_state(npf, elem, node, th, ud):
    """
    Compute hydrostatic background state for atmospheric model.
    Handles arbitrary stratification profiles and proper numerical integration.
    """
    # Thermodynamic constants
    Gamma = th.gm1 / th.gamm
    gamm = th.gamm
    gm1 = th.gm1
    Gamma_inv = 1.0 / Gamma
    gm1_inv = 1.0 / gm1

    # Grid parameters along the vertical axis
    vv = axes.vertical_axis(ud)
    icy = elem.sc[vv]
    igy = elem.igs[vv]
    dyv = elem.dxyz[vv]
    y_c = axes.coords_along(elem, vv)
    y_n = axes.coords_along(node, vv)

    # Reference state at y=0
    rhoY0 = 1.0
    g = ud.gravity_strength[axes.vertical_axis(ud)]
    p0 = rhoY0**gamm
    pi0 = rhoY0**gm1

    if g != 0.0 and elem.metric is not None:
        _integrated_state_fields(npf, elem, node, th, ud)
        return

    if g != 0.0:
        ###########################
        # Update cell hydrostates
        ###########################

        # Define midpoint quadrature along vertical (y-axis)
        dys = np.hstack(
            (
                np.ones(igy - 1) * -dyv,
                [-dyv / 2],
                [dyv / 2],
                np.ones(icy - 3) * dyv,
            )
        )

        # Cell centers and midpoints for integration
        y_ps = y_c
        y_ms = np.hstack((y_c[1:igy], y_n[igy], y_n[igy], y_c[igy:-1]))

        # Get inverse stratification at each point
        S_ps = 1.0 / ud.stratification(y_ps)
        S_ms = 1.0 / ud.stratification(y_ms)

        # Trapezoidal integration over inverse stratification
        S_integral_p = 0.5 * dys * (S_ms + S_ps)

        # Cumulative integration (split at boundary igy)
        S_integral_p[:igy] = np.cumsum(S_integral_p[:igy][::-1])[::-1]
        S_integral_p[igy:] = np.cumsum(S_integral_p[igy:])

        # Calculate hydrostatic fields
        pi_hydro = pi0 - Gamma * g * S_integral_p
        p_hydro = pi_hydro**Gamma_inv
        rhoY_hydro = pi_hydro**gm1_inv

        # Update cell solutions
        npf.HydroState.rhoY0[:] = rhoY_hydro
        npf.HydroState.rho0[:] = rhoY_hydro * S_ps
        npf.HydroState.p0[:] = p_hydro
        npf.HydroState.p20[:] = pi_hydro / ud.Msq
        npf.HydroState.S0[:] = S_ps
        npf.HydroState.S10[:] = 0.0
        npf.HydroState.Y0[:] = 1.0 / S_ps

        ############################
        # Update node hydrostates
        ############################

        # Bottom reference node (y=0)
        npf.HydroState_n.Y0[igy] = ud.stratification(0.0)
        npf.HydroState_n.rhoY0[igy] = rhoY0
        npf.HydroState_n.rho0[igy] = rhoY0 / ud.stratification(0.0)
        npf.HydroState_n.S0[igy] = 1.0 / npf.HydroState_n.Y0[igy]
        npf.HydroState_n.p0[igy] = p0
        npf.HydroState_n.p20[igy] = pi0 / ud.Msq

        # Ghost cells below bottom (negative heights)
        Sn_integral_p = np.zeros(igy)
        yn_p = y_n[:igy] - dyv
        yn_m = y_n[1 : igy + 1] - dyv

        Sn_integral_p[:] = -dyv * 1.0 / ud.stratification(0.5 * (yn_p + yn_m))
        Sn_integral_p = np.cumsum(Sn_integral_p[:igy][::-1])[::-1]

        # Bulk domain above reference level
        yn_p = y_n[igy + 1 :]
        yn_m = np.zeros_like(yn_p)
        yn_m[1:] = yn_p[:-1]

        Sn_p = 1.0 / ud.stratification(0.5 * (yn_p + yn_m))
        Sn_integral_p = np.hstack((Sn_integral_p, np.cumsum(dyv * Sn_p)))

        # Calculate nodal hydrostatic fields
        pi_hydro_n = pi0 - Gamma * g * Sn_integral_p
        rhoY_hydro_n = pi_hydro_n**gm1_inv

        # Update node solutions - below reference
        npf.HydroState_n.rhoY0[:igy] = rhoY_hydro_n[:igy]
        npf.HydroState_n.Y0[: igy + 1] = ud.stratification(
            0.5 * (y_ps[: igy + 1] + y_ps[: igy + 1] - dyv)
        )
        npf.HydroState_n.rho0[:igy] = rhoY_hydro_n[:igy] / npf.HydroState_n.Y0[:igy]
        npf.HydroState_n.S0[:igy] = 1.0 / npf.HydroState_n.Y0[:igy]
        npf.HydroState_n.p0[:igy] = rhoY_hydro_n[:igy] ** th.gamm
        npf.HydroState_n.p20[:igy] = pi_hydro_n[:igy] / ud.Msq

        # Update node solutions - above reference
        npf.HydroState_n.rhoY0[igy + 1 :] = rhoY_hydro_n[igy:]
        npf.HydroState_n.Y0[igy + 1 :] = ud.stratification(
            0.5 * (y_ps[igy:] + y_ps[igy:] + dyv)
        )
        npf.HydroState_n.rho0[igy + 1 :] = (
            rhoY_hydro_n[igy:] / npf.HydroState_n.Y0[igy + 1 :]
        )
        npf.HydroState_n.S0[igy + 1 :] = 1.0 / npf.HydroState_n.Y0[igy + 1 :]
        npf.HydroState_n.p0[igy + 1 :] = rhoY_hydro_n[igy:] ** th.gamm
        npf.HydroState_n.p20[igy + 1 :] = pi_hydro_n[igy:] / ud.Msq

    else:
        # No gravity case - uniform atmosphere
        npf.HydroState.p20[:] = 1.0
        npf.HydroState.p0[:] = 1.0
        npf.HydroState.rho0[:] = 1.0
        npf.HydroState.rhoY0[:] = 1.0
        npf.HydroState.Y0[:] = 1.0
        npf.HydroState.S0[:] = 1.0
        npf.HydroState.S10[:] = 0.0

        npf.HydroState_n.p20[:] = 1.0
        npf.HydroState_n.p0[:] = 1.0
        npf.HydroState_n.rho0[:] = 1.0
        npf.HydroState_n.rhoY0[:] = 1.0
        npf.HydroState_n.Y0[:] = 1.0
        npf.HydroState_n.S0[:] = 1.0


def _integrated_state_fields(npf, elem, node, th, ud):
    """Terrain branch of integrated_state: hydrostates as per-column fields.

    The Exner pressure follows from quadrature of the inverse stratification
    on a fine auxiliary 1D z-grid (the background depends on physical height
    only), evaluated at the cell/node height fields z(xi, eta) by linear
    interpolation. The pi reference (rhoY = 1) sits at z = 0, matching the
    profile branch.
    """
    Gamma = th.gm1 / th.gamm
    Gamma_inv = 1.0 / Gamma
    gm1_inv = 1.0 / th.gm1

    vv = axes.vertical_axis(ud)
    g = ud.gravity_strength[vv]

    # generalized altitude: == z for vertical-line maps, r - a on a sphere
    z_c = elem.metric.height
    z_n = node.metric.height
    z_lo = min(z_c.min(), z_n.min(), 0.0)
    z_hi = max(z_c.max(), z_n.max(), 0.0)
    nfine = max(_HYDRO_QUAD_MIN_POINTS, _HYDRO_QUAD_POINTS_PER_CELL * int(elem.sc[vv]))
    zf = np.linspace(z_lo, z_hi, nfine)

    Sf = 1.0 / ud.stratification(zf)
    integral = np.concatenate(
        ([0.0], np.cumsum(0.5 * (Sf[1:] + Sf[:-1]) * np.diff(zf)))
    )
    integral -= np.interp(0.0, zf, integral)

    rhoY0 = 1.0
    pi0 = rhoY0**th.gm1

    for states, z in ((npf.HydroState, z_c), (npf.HydroState_n, z_n)):
        pi = pi0 - Gamma * g * np.interp(z, zf, integral)
        S = 1.0 / ud.stratification(z)
        states.rhoY0[...] = pi**gm1_inv
        states.p0[...] = pi**Gamma_inv
        states.p20[...] = pi / ud.Msq
        states.S0[...] = S
        states.S10[...] = 0.0
        states.Y0[...] = 1.0 / S
        states.rho0[...] = states.rhoY0 * S


def analytical_state(npf, elem, node, th, ud):
    """Isothermal hydrostatic background, discrete-exact per cell.

    With terrain the same closed form is evaluated at the physical heights
    z(xi, eta) with the local vertical cell extent dz = J * deta, so the
    hydrostates become full per-column fields (States in field mode);
    without terrain the expressions reduce to the legacy 1D profiles
    bit-identically.
    """
    vv = axes.vertical_axis(ud)
    g = ud.gravity_strength[vv]
    Gamma = th.Gamma
    Hex = 1.0 / (th.Gamma * g)
    dy = elem.dxyz[vv]

    if elem.metric is not None and not elem.metric.vertical_line:
        # general map (sphere): altitude is metric.height (r - a) and the
        # vertical arc length per unit eta is |t_v| = metric.h_v
        mn, mc = node.metric, elem.metric
        z_n, dz_n = mn.height, mn.h_v * dy
        z_c, dz_c = mc.height, mc.h_v * dy
    elif elem.metric is not None:
        # local vertical cell extent dz = z_eta * deta with z_eta = J/(N_v)_v:
        # on horizontally stretched grids J = x' z_eta y' is the VOLUME
        # measure, not the height increment; for vertical-line maps
        # (N_v)_v == 1 and this is bit-exactly the legacy J * deta
        mn, mc = node.metric, elem.metric
        z_eta_n = mn.J / mn.N[mn.vaxis][mn.cart_v]
        z_eta_c = mc.J / mc.N[mc.vaxis][mc.cart_v]
        z_n, dz_n = mn.z, z_eta_n * dy
        z_c, dz_c = mc.z, z_eta_c * dy
    else:
        z_n, dz_n = axes.coords_along(node, vv), dy
        z_c, dz_c = axes.coords_along(elem, vv), dy

    pi_np = np.exp(-(z_n + 0.5 * dz_n) / Hex)
    pi_nm = np.exp(-(z_n - 0.5 * dz_n) / Hex)
    pi_n = np.exp(-(z_n) / Hex)

    Y_n = -Gamma * g * dz_n / (pi_np - pi_nm)
    P_n = pi_n**th.gm1inv
    p_n = pi_n**th.Gammainv
    rho_n = P_n / Y_n

    npf.HydroState_n.p20[...] = pi_n / ud.Msq
    npf.HydroState_n.p0[...] = p_n
    npf.HydroState_n.rho0[...] = rho_n
    npf.HydroState_n.rhoY0[...] = P_n
    npf.HydroState_n.Y0[...] = Y_n
    npf.HydroState_n.S0[...] = 1.0 / Y_n

    pi_cp = np.exp(-(z_c + 0.5 * dz_c) / Hex)
    pi_cm = np.exp(-(z_c - 0.5 * dz_c) / Hex)
    pi_c = np.exp(-(z_c) / Hex)

    Y_c = -Gamma * g * dz_c / (pi_cp - pi_cm)
    P_c = pi_c**th.gm1inv
    p_c = pi_c**th.Gammainv
    rho_c = P_c / Y_c

    npf.HydroState.p20[...] = pi_c / ud.Msq
    npf.HydroState.p0[...] = p_c
    npf.HydroState.rho0[...] = rho_c
    npf.HydroState.rhoY0[...] = P_c
    npf.HydroState.Y0[...] = Y_c
    npf.HydroState.S0[...] = 1.0 / Y_c


def initial_pressure(Sol, npf, elem, node, ud, th):
    """2D x-y initial-condition helper (vertical = axis 1 by convention)."""
    assert elem.ndim == 2, "initial_pressure() is a 2D x-y IC helper"
    assert elem.metric is None, "initial_pressure() does not support terrain"
    Gammainv = th.Gammainv
    igy = node.igy
    igx = node.igx
    icx = elem.icx
    dx = node.dx
    dy = node.dy

    beta = np.zeros((node.icx))
    bdpdx = np.zeros((node.icx))

    x_idx_m = slice(0, -1)
    x_idx_c = slice(1, None)
    y_idx = slice(igy, -igy)
    xn_idx = slice(1, -1)
    height = node.y[-igy - 1]

    Pc = Sol.rhoY[x_idx_c, y_idx]
    Pm = Sol.rhoY[x_idx_m, y_idx]
    thc = Pc / Sol.rho[x_idx_c, y_idx]
    thm = Pm / Sol.rho[x_idx_m, y_idx]
    beta[xn_idx] = np.sum(0.5 * (Pm * thm + Pc * thc) * dy, axis=1)
    bdpdx[xn_idx] = np.sum(
        0.5
        * (Pm * thm + Pc * thc)
        * (npf.p2_cells[x_idx_c, y_idx] - npf.p2_cells[x_idx_m, y_idx])
        * dy,
        axis=1,
    )

    beta *= Gammainv / height
    bdpdx *= Gammainv / height / dx

    coeff = np.zeros((elem.icx))
    pibot = np.zeros((elem.icx))
    coeff[igx + 1 : -igx + 1] = np.cumsum(coeff[igx:-igx] + dx / beta[igx + 1 : -igx])
    pibot[igx + 1 : -igx + 1] = np.cumsum(
        pibot[igx:-igx] - dx * bdpdx[igx + 1 : -igx] / beta[igx + 1 : -igx]
    )

    dotPU = pibot[icx - igx] / coeff[icx - igx]
    pibot[igx:-igx] -= dotPU * coeff[igx:-igx]

    x_idx = slice(igx, -igx + 1)
    y_idx = slice(igy, -igy + 1)

    npf.p2_cells[x_idx, y_idx] += pibot[x_idx].reshape(
        -1, 1
    ) - 1.0 * npf.HydroState.p20[y_idx].reshape(1, -1)

    icxn = node.icx
    icyn = node.icy
    x_idx = slice(1, icxn - 1)
    y_idx = slice(igy, -igy + 1)
    height = node.y[-igy]

    Pc = Sol.rhoY[1:, y_idx]
    thc = Pc / Sol.rho[1:, y_idx]

    beta = np.zeros((elem.icx,))
    bdpdx = np.zeros((elem.icx))

    beta[1:] = np.sum(Pc * thc * dy, axis=1)
    beta *= Gammainv / height

    bdpdx[1:] = np.sum(
        Pc * thc * (npf.p2_nodes[1:-1, igy:-igy] - npf.p2_nodes[:-2, igy:-igy]) * dy,
        axis=1,
    )
    bdpdx *= Gammainv / height / dx

    coeff = np.zeros((node.icx))
    pibot = np.zeros((node.icx))

    coeff[igx + 1 : -igx + 1] = np.cumsum(coeff[igx:-igx] + dx / beta[igx + 1 :])
    pibot[igx + 1 : -igx + 1] = np.cumsum(
        pibot[igx:-igx] - dx * bdpdx[igx + 1 :] / beta[igx + 1 :]
    )

    dotPU = pibot[icx - igx] / coeff[icx - igx]

    pibot[igx:-igx] -= dotPU * coeff[igx:-igx]

    x_idx = slice(igx, -igx + 1)
    y_idx = slice(igy, -igy + 1)
    npf.p2_nodes[x_idx, y_idx] += pibot[x_idx].reshape(
        -1, 1
    ) - 1.0 * npf.HydroState_n.p20[y_idx].reshape(1, -1)

    npf.dp2_nodes[:, :] = npf.p2_nodes

    # guess initial node value (at left-most node)
    npf.p2_nodes[igx, igy:-igy] = npf.dp2_nodes[igx, igy:-igy]

    npf.p2_nodes[:, :] = __loop_over_array(
        igx, igy, icxn, icyn, npf.p2_nodes, npf.dp2_nodes
    )

    assert ((node.icx + 1) % 2) == 1
    delp2 = 0.5 * (npf.p2_nodes[-igx - 1, igy:-igy] - npf.p2_nodes[igx, igy:-igy])
    delp2 = delp2.reshape(1, -1)
    sgn = np.ones_like(npf.p2_nodes[:, 0][igy:-igy]).reshape(-1, 1)

    sgn[1::2] *= -1

    npf.p2_nodes[igx:-igx, igy:-igy] += sgn * delp2
    bdry_n.set_ghost_nodes(npf.p2_nodes, node, ud)

    npf.dp2_nodes[:, :] = 0.0

    inner_domain = (slice(igx, -igx), slice(igy, -igy))
    pi = ud.Msq * (npf.p2_cells[inner_domain] + 1.0 * npf.HydroState.p20[igy:-igy])
    Y = Sol.rhoY[inner_domain] / Sol.rho[inner_domain]
    rhoold = np.copy(Sol.rho[inner_domain])
    Sol.rhoY[inner_domain] = pi**th.gm1inv
    Sol.rho[inner_domain] = Sol.rhoY[inner_domain] / Y
    Sol.rhou[inner_domain] *= Sol.rho[inner_domain] / rhoold
    Sol.rhov[inner_domain] *= Sol.rho[inner_domain] / rhoold
    Sol.rhow[inner_domain] *= Sol.rho[inner_domain] / rhoold
    Sol.rhoX[inner_domain] *= Sol.rho[inner_domain] / rhoold


# need details:
# populate the rest of the nodes recursively based on the left-most node.
# recursive: use numba.
@nb.jit(nopython=True)
def __loop_over_array(igx, igy, icxn, icyn, p, dp):
    for j in range(igy, icyn - igy):
        for i in range(igx + 1, icxn - igx):
            p[i, j] = 2.0 * dp[i - 1, j] - p[i - 1, j]
    return p
