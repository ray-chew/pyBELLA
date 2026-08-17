import numpy as np
import numba as nb

from ..utils import options as opts

from ..flow_solver.utils import fields
from ..flow_solver.utils.boundary import node_boundary as bdry_n
from ..flow_solver.physics import hydrostatics

from .case_setup import build_bdry, make_diag_state


class UserData(object):
    # planetary -> 160.0;  long-wave -> 20.0;  standard -> 1.0;
    scale_factor = 20.0

    def __init__(self):
        self.scale_factor = self.scale_factor

        self.h_ref = 10000.0  # [m]
        self.t_ref = 100.0  # [s]
        self.T_ref = 300.00  # [K]
        self.p_ref = 1e5  # [Pa]
        self.omega = 7.292 * 1e-5  # [s^{-1}]
        self.grav = 9.81  # [m/s^2]
        self.R_gas = 287.4  # [J kg^{-1} K^{-1}]
        self.u_ref = self.h_ref / self.t_ref  # [m/s]
        self.Nsq_ref = 1.0e-4  # [s^{-2}]
        self.Msq = (
            self.u_ref * self.u_ref / (self.R_gas * self.T_ref)
        )  # Mach number squared

        self.gravity_strength = np.zeros((3))

        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)

        gravity_mask = (self.gravity_strength > np.finfo(np.float64).eps) | (
            np.arange(3) == 1
        )
        self.i_gravity = gravity_mask.astype(int)
        if np.any(gravity_mask):
            self.gravity_direction = np.where(gravity_mask)[0][
                -1
            ]  # Use last matching index

        self.xmin = -15.0 * self.scale_factor
        self.xmax = 15.0 * self.scale_factor
        self.ymin = 0.0
        self.ymax = 1.0
        self.zmin = -1.0
        self.zmax = 1.0

        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.PERIODIC
        )

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9

        self.dtfixed0 = (
            10.0 * (12.5 / 15.0) * 0.5 * self.scale_factor * 30.0 / self.t_ref
        )
        self.dtfixed = (
            10.0 * (12.5 / 15.0) * 0.5 * self.scale_factor * 30.0 / self.t_ref
        )
        # self.dtfixed0 = 5.0 * (12.5 / 15.0) * 0.5 * self.scale_factor * 30.0 / self.t_ref
        # self.dtfixed = 5.0 * (12.5 / 15.0) * 0.5 * self.scale_factor * 30.0 / self.t_ref

        self.dtfixed0 = 1.0
        self.dtfixed = 1.0

        self.inx = 301 + 1
        self.iny = 10 + 1
        self.inz = 1

        self.tout = [self.scale_factor * 1.0 * 3000.0 / self.t_ref]

        self.tol = 1.0e-12
        self.stepmax = 31
        self.max_iterations = 6000

        self.autogen_fn = False

        self.output_timesteps = True

        self.stratification = self.stratification_function
        self.molly = self.molly_function
        self.rhoe = self.rhoe_method

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_internal_long_wave"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.diag_state = make_diag_state(
            "test_internal_long_wave",
            "target_internal_long_wave",
            self.inx,
            self.iny,
            self.stepmax,
        )

    def stratification_function(self, y):
        Nsq = self.Nsq_ref * self.t_ref * self.t_ref
        g = self.gravity_strength[1] / self.Msq

        return np.exp(Nsq * y / g)

    def molly_function(self, x):
        del0 = 0.25
        L = self.xmax - self.xmin
        xi_l = np.minimum(1.0, (x - self.xmin) / (del0 * L))
        xi_r = np.minimum(1.0, (self.xmax - x) / (del0 * L))

        return 0.5 * np.minimum(1.0 - np.cos(np.pi * xi_l), 1.0 - np.cos(np.pi * xi_r))

    @staticmethod
    def rhoe_method(rho, u, v, w, p, ud, th):
        Msq = ud.compressibility * ud.Msq

        gm1inv = th.gm1inv
        return p * gm1inv + 0.5 * Msq * rho * (u * u + v * v + w * w)


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed
    delth = 0.01 / ud.T_ref
    xc = -0.0 * ud.scale_factor * 50.0e3 / ud.h_ref
    xc = 0.0
    a = ud.scale_factor * 5.0e3 / ud.h_ref

    # ensemble-member perturbation of the theta' wave (3D DA OSSE):
    # seeded amplitude/position shifts + a z-modulation used by the 3D
    # branch; inert for every shipped config (seed None / perturb_type
    # default 'pos_perturb')
    zamp, zphase = 0.0, 0.0
    if seed is not None and ud.perturb_type == "igw_theta":
        np.random.seed(seed)
        xi = np.random.random(3)
        delth *= 1.0 + 0.4 * (xi[0] - 0.5)
        xc += 2.0 * a * (xi[1] - 0.5)
        zamp, zphase = 0.2, 2.0 * np.pi * xi[2]

    hydrostatics.integrated_state(npf, elem, node, th, ud)

    if elem.ndim == 3:
        return _sol_init_3d(
            Sol, npf, elem, node, th, ud, u0, v0, w0, delth, xc, a, zamp, zphase
        )

    HySt = fields.States(node.sc)
    HyStn = fields.States(node.sc)

    x = elem.x.reshape(-1, 1)
    y = elem.y.reshape(1, -1)

    Y = ud.stratification(y) + delth * ud.molly(x) * np.sin(np.pi * y) / (
        1.0 + (x - xc) ** 2 / (a**2)
    )

    xn = node.x[:-1].reshape(-1, 1)
    yn = node.y[:-1].reshape(1, -1)

    Yn = ud.stratification(yn) + delth * ud.molly(xn) * np.sin(np.pi * yn) / (
        1.0 + (xn - xc) ** 2 / (a**2)
    )

    hydrostatics.column(HySt, HyStn, Y, Yn, elem, node, th, ud)

    x_idx = slice(None)
    y_idx = slice(elem.igy, -elem.igy + 1)
    xc_idx = slice(0, -1)
    yc_idx = slice(0, -1)
    c_idx = (xc_idx, yc_idx)

    u, v, w = u0, v0, w0
    if ud.is_compressible:
        p = HySt.p0[:, y_idx][c_idx]
        rhoY = HySt.rhoY0[:, y_idx][c_idx]
    else:
        p = npf.HydroState.p0[y_idx]
        rhoY = npf.HydroState.rhoY0[y_idx]

    rho = rhoY / Y[:, y_idx]
    Sol.rho[x_idx, y_idx] = rho
    Sol.rhou[x_idx, y_idx] = rho * u
    Sol.rhov[x_idx, y_idx] = rho * v
    Sol.rhow[x_idx, y_idx] = rho * w
    Sol.rhoY[x_idx, y_idx] = rhoY

    npf.p2_cells[x_idx, y_idx] = HySt.p20[x_idx, y_idx][c_idx]

    Sol.rhoX[x_idx, y_idx] = Sol.rho[x_idx, y_idx] * (
        1.0 / Y[:, y_idx] - npf.HydroState.S0[y_idx]
    )

    npf.p2_nodes[:, elem.igy : -elem.igy] = HyStn.p20[:, elem.igy : -elem.igy]

    hydrostatics.initial_pressure(Sol, npf, elem, node, ud, th)

    ud.nonhydrostasy = 1.0 if ud.is_nonhydrostatic == 1 else 0.0
    ud.compressibility = 1.0 if ud.is_compressible == 1 else 0.0

    if "imbal" in ud.aux:
        npf.p2_nodes[...] = 0.0

    return Sol


# ---------------------------------------------------------------------------
# 3D branch (used by the 3D data-assimilation OSSE). Faithful transcriptions
# of the 2D path with z as a trailing broadcast axis: with zamp == 0 every
# z-slab of every field is BITWISE identical to the 2D IC
# (test_scripts/test_igw3d_ic_oracle.py).
# hydrostatics.column/initial_pressure stay untouched (they assert ndim == 2);
# the twins live here because they are igw-IC helpers, not solver code.
# ---------------------------------------------------------------------------


def _sol_init_3d(Sol, npf, elem, node, th, ud, u0, v0, w0, delth, xc, a, zamp, zphase):
    HySt = fields.States(node.sc)
    HyStn = fields.States(node.sc)

    x = elem.x.reshape(-1, 1, 1)
    y = elem.y.reshape(1, -1, 1)
    z = elem.z.reshape(1, 1, -1)

    xn = node.x[:-1].reshape(-1, 1, 1)
    yn = node.y[:-1].reshape(1, -1, 1)
    zn = node.z[:-1].reshape(1, 1, -1)

    Lz = ud.zmax - ud.zmin
    zmod = 1.0 + zamp * np.cos(2.0 * np.pi * (z - ud.zmin) / Lz + zphase)
    zmod_n = 1.0 + zamp * np.cos(2.0 * np.pi * (zn - ud.zmin) / Lz + zphase)

    Y = ud.stratification(y) + delth * zmod * ud.molly(x) * np.sin(np.pi * y) / (
        1.0 + (x - xc) ** 2 / (a**2)
    )
    Yn = ud.stratification(yn) + delth * zmod_n * ud.molly(xn) * np.sin(np.pi * yn) / (
        1.0 + (xn - xc) ** 2 / (a**2)
    )

    _column3(HySt, HyStn, Y, Yn, elem, node, th, ud)

    x_idx = slice(None)
    y_idx = slice(elem.igy, -elem.igy + 1)
    z_idx = slice(None)
    c_idx = (slice(0, -1), slice(0, -1), slice(0, -1))

    u, v, w = u0, v0, w0
    if ud.is_compressible:
        p = HySt.p0[:, y_idx, :][c_idx]
        rhoY = HySt.rhoY0[:, y_idx, :][c_idx]
    else:
        p = npf.HydroState.p0[y_idx].reshape(1, -1, 1)
        rhoY = npf.HydroState.rhoY0[y_idx].reshape(1, -1, 1)

    rho = rhoY / Y[:, y_idx, :]
    Sol.rho[x_idx, y_idx, z_idx] = rho
    Sol.rhou[x_idx, y_idx, z_idx] = rho * u
    Sol.rhov[x_idx, y_idx, z_idx] = rho * v
    Sol.rhow[x_idx, y_idx, z_idx] = rho * w
    Sol.rhoY[x_idx, y_idx, z_idx] = rhoY

    npf.p2_cells[x_idx, y_idx, z_idx] = HySt.p20[x_idx, y_idx, z_idx][c_idx]

    Sol.rhoX[x_idx, y_idx, z_idx] = Sol.rho[x_idx, y_idx, z_idx] * (
        1.0 / Y[:, y_idx, :] - npf.HydroState.S0[y_idx].reshape(1, -1, 1)
    )

    npf.p2_nodes[:, elem.igy : -elem.igy, :] = HyStn.p20[:, elem.igy : -elem.igy, :]

    _initial_pressure3(Sol, npf, elem, node, ud, th)

    ud.nonhydrostasy = 1.0 if ud.is_nonhydrostatic == 1 else 0.0
    ud.compressibility = 1.0 if ud.is_compressible == 1 else 0.0

    if "imbal" in ud.aux:
        npf.p2_nodes[...] = 0.0

    return Sol


def _column3(HydroState, HydroState_n, Y, Y_n, elem, node, th, ud):
    """hydrostatics.column twin: vertical ops along axis 1, trailing z."""
    Gamma = th.gm1 / th.gamm
    gamm = th.gamm
    gm1 = th.gm1
    Gamma_inv = 1.0 / Gamma
    gm1_inv = 1.0 / gm1

    icy = elem.icy
    igy = elem.igy

    c_idx = (slice(0, -1), slice(0, -1), slice(0, -1))
    xzc_idx = (slice(0, -1), slice(0, -1))

    rhoY0 = 1.0

    from ..utils import axes

    g = ud.gravity_strength[axes.vertical_axis(ud)]

    p0 = rhoY0**gamm
    pi0 = rhoY0**gm1
    xc_idx = slice(0, -1)
    zc_idx = slice(0, -1)
    HydroState_n.rho0[xc_idx, igy, zc_idx] = rhoY0 / Y_n[:, igy, :]
    HydroState_n.rhoY0[xc_idx, igy, zc_idx] = rhoY0
    HydroState_n.Y0[xc_idx, igy, zc_idx] = Y_n[:, igy, :]
    HydroState_n.S0[xc_idx, igy, zc_idx] = 1.0 / Y_n[:, igy, :]
    HydroState_n.p0[xc_idx, igy, zc_idx] = p0
    HydroState_n.p20[xc_idx, igy, zc_idx] = pi0 / ud.Msq

    dys = np.array(
        [-elem.dy] + [-elem.dy / 2] + [elem.dy / 2] + list(np.ones((icy - 3)) * elem.dy)
    ).reshape(1, -1, 1)
    S_p = 1.0 / Y[:, :, :]
    S_m = np.zeros_like(S_p)
    S_m[:, igy - 1 : igy + 1, :] = (1.0 / Y_n[:, igy, :])[:, None, :]
    S_m[:, 0, :] = 1.0 / Y[:, igy - 1, :]
    S_m[:, igy + 1 :, :] = 1.0 / Y[:, igy:-1, :]

    S_integral_p = dys * 0.5 * (S_p + S_m)
    S_integral_p[:, :igy, :] = np.cumsum(S_integral_p[:, :igy, :][:, ::-1, :], axis=1)[
        :, ::-1, :
    ]
    S_integral_p[:, igy:, :] = np.cumsum(S_integral_p[:, igy:, :], axis=1)

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

    Sn_p = 1.0 / Y[:, :, :]
    dys = (np.ones((icy)) * elem.dy).copy()
    dys[:igy] *= -1
    dys = dys.reshape(1, -1, 1)
    Sn_integral_p = dys * Sn_p
    Sn_integral_p[:, :igy, :] = np.cumsum(
        Sn_integral_p[:, :igy, :][:, ::-1, :], axis=1
    )[:, ::-1, :]
    Sn_integral_p[:, igy:, :] = np.cumsum(Sn_integral_p[:, igy:, :], axis=1)

    pi_hydro_n = pi0 - Gamma * g * Sn_integral_p
    rhoY_hydro_n = pi_hydro_n**gm1_inv

    HydroState_n.rhoY0[xc_idx, :igy, zc_idx] = rhoY_hydro_n[:, :igy, :]
    HydroState_n.Y0[xc_idx, :igy, zc_idx] = Y_n[0, :igy, :]
    HydroState_n.S0[xc_idx, :igy, zc_idx] = 1.0 / Y_n[:, :igy, :]
    HydroState_n.p0[xc_idx, :igy, zc_idx] = rhoY_hydro_n[:, :igy, :] ** th.gamm
    HydroState_n.p20[xc_idx, :igy, zc_idx] = pi_hydro_n[:, :igy, :] / ud.Msq

    HydroState_n.rhoY0[xc_idx, igy + 1 :, zc_idx] = rhoY_hydro_n[:, igy:, :]
    HydroState_n.Y0[xc_idx, igy + 1 :, zc_idx] = Y_n[0, igy:, :]
    HydroState_n.S0[xc_idx, igy + 1 :, zc_idx] = 1.0 / Y_n[:, igy:, :]
    HydroState_n.p0[xc_idx, igy + 1 :, zc_idx] = rhoY_hydro_n[:, igy:, :] ** th.gamm
    HydroState_n.p20[xc_idx, igy + 1 :, zc_idx] = pi_hydro_n[:, igy:, :] / ud.Msq


def _initial_pressure3(Sol, npf, elem, node, ud, th):
    """hydrostatics.initial_pressure twin: the x-line pressure corrections run
    per z-plane (trailing axis). z-pairing convention: line integrals are
    built from CELL fields on cell z-planes; node arrays enter/receive on
    their LOWER-node z-plane ([:, :, :-1]) — exact for z-uniform slabs (the
    IC oracle checks bitwise), an O(dz^2)-consistent per-plane approximation
    for the weakly z-modulated members."""
    Gammainv = th.Gammainv
    igy = node.igy
    igx = node.igx
    icx = elem.icx
    icz = elem.icz
    dx = node.dx
    dy = node.dy

    beta = np.zeros((node.icx, icz))
    bdpdx = np.zeros((node.icx, icz))

    x_idx_m = slice(0, -1)
    x_idx_c = slice(1, None)
    y_idx = slice(igy, -igy)
    xn_idx = slice(1, -1)
    height = node.y[-igy - 1]

    Pc = Sol.rhoY[x_idx_c, y_idx, :]
    Pm = Sol.rhoY[x_idx_m, y_idx, :]
    thc = Pc / Sol.rho[x_idx_c, y_idx, :]
    thm = Pm / Sol.rho[x_idx_m, y_idx, :]
    beta[xn_idx] = np.sum(0.5 * (Pm * thm + Pc * thc) * dy, axis=1)
    bdpdx[xn_idx] = np.sum(
        0.5
        * (Pm * thm + Pc * thc)
        * (npf.p2_cells[x_idx_c, y_idx, :] - npf.p2_cells[x_idx_m, y_idx, :])
        * dy,
        axis=1,
    )

    beta *= Gammainv / height
    bdpdx *= Gammainv / height / dx

    coeff = np.zeros((elem.icx, icz))
    pibot = np.zeros((elem.icx, icz))
    coeff[igx + 1 : -igx + 1] = np.cumsum(
        coeff[igx:-igx] + dx / beta[igx + 1 : -igx], axis=0
    )
    pibot[igx + 1 : -igx + 1] = np.cumsum(
        pibot[igx:-igx] - dx * bdpdx[igx + 1 : -igx] / beta[igx + 1 : -igx], axis=0
    )

    dotPU = pibot[icx - igx] / coeff[icx - igx]
    pibot[igx:-igx] -= dotPU * coeff[igx:-igx]

    x_idx = slice(igx, -igx + 1)
    y_idx = slice(igy, -igy + 1)

    npf.p2_cells[x_idx, y_idx, :] += pibot[x_idx][
        :, None, :
    ] - 1.0 * npf.HydroState.p20[y_idx].reshape(1, -1, 1)

    icxn = node.icx
    icyn = node.icy
    iczn = node.icz
    x_idx = slice(1, icxn - 1)
    y_idx = slice(igy, -igy + 1)
    height = node.y[-igy]

    Pc = Sol.rhoY[1:, y_idx, :]
    thc = Pc / Sol.rho[1:, y_idx, :]

    beta = np.zeros((elem.icx, icz))
    bdpdx = np.zeros((elem.icx, icz))

    beta[1:] = np.sum(Pc * thc * dy, axis=1)
    beta *= Gammainv / height

    bdpdx[1:] = np.sum(
        Pc
        * thc
        * (npf.p2_nodes[1:-1, igy:-igy, :-1] - npf.p2_nodes[:-2, igy:-igy, :-1])
        * dy,
        axis=1,
    )
    bdpdx *= Gammainv / height / dx

    coeff = np.zeros((node.icx, icz))
    pibot = np.zeros((node.icx, icz))

    coeff[igx + 1 : -igx + 1] = np.cumsum(
        coeff[igx:-igx] + dx / beta[igx + 1 :], axis=0
    )
    pibot[igx + 1 : -igx + 1] = np.cumsum(
        pibot[igx:-igx] - dx * bdpdx[igx + 1 :] / beta[igx + 1 :], axis=0
    )

    dotPU = pibot[icx - igx] / coeff[icx - igx]

    pibot[igx:-igx] -= dotPU * coeff[igx:-igx]

    # node arrays receive the correction on their lower-node z-plane; the top
    # node plane (periodic z) mirrors plane 0
    pibot_n = np.concatenate([pibot, pibot[:, :1]], axis=1)

    x_idx = slice(igx, -igx + 1)
    y_idx = slice(igy, -igy + 1)
    npf.p2_nodes[x_idx, y_idx, :] += pibot_n[x_idx][
        :, None, :
    ] - 1.0 * npf.HydroState_n.p20[y_idx].reshape(1, -1, 1)

    npf.dp2_nodes[...] = npf.p2_nodes

    # guess initial node value (at left-most node)
    npf.p2_nodes[igx, igy:-igy, :] = npf.dp2_nodes[igx, igy:-igy, :]

    npf.p2_nodes[...] = _loop_over_array3(
        igx, igy, icxn, icyn, iczn, npf.p2_nodes, npf.dp2_nodes
    )

    assert ((node.icx + 1) % 2) == 1
    delp2 = 0.5 * (npf.p2_nodes[-igx - 1, igy:-igy, :] - npf.p2_nodes[igx, igy:-igy, :])
    delp2 = delp2[None, :, :]
    sgn = np.ones((npf.p2_nodes.shape[0] - 2 * igx, 1, 1))
    sgn[1::2] *= -1

    npf.p2_nodes[igx:-igx, igy:-igy, :] += sgn * delp2
    bdry_n.set_ghost_nodes(npf.p2_nodes, node, ud)

    npf.dp2_nodes[...] = 0.0

    inner_domain = (slice(igx, -igx), slice(igy, -igy), slice(None))
    pi = ud.Msq * (
        npf.p2_cells[inner_domain]
        + 1.0 * npf.HydroState.p20[igy:-igy].reshape(1, -1, 1)
    )
    Y = Sol.rhoY[inner_domain] / Sol.rho[inner_domain]
    rhoold = np.copy(Sol.rho[inner_domain])
    Sol.rhoY[inner_domain] = pi**th.gm1inv
    Sol.rho[inner_domain] = Sol.rhoY[inner_domain] / Y
    Sol.rhou[inner_domain] *= Sol.rho[inner_domain] / rhoold
    Sol.rhov[inner_domain] *= Sol.rho[inner_domain] / rhoold
    Sol.rhow[inner_domain] *= Sol.rho[inner_domain] / rhoold
    Sol.rhoX[inner_domain] *= Sol.rho[inner_domain] / rhoold


@nb.jit(nopython=True)
def _loop_over_array3(igx, igy, icxn, icyn, iczn, p, dp):
    for k in range(iczn):
        for j in range(igy, icyn - igy):
            for i in range(igx + 1, icxn - igx):
                p[i, j, k] = 2.0 * dp[i - 1, j, k] - p[i - 1, j, k]
    return p
