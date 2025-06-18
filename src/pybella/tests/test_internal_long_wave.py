import numpy as np

from ..utils import options as opts

from ..flow_solver.utils import fields
from ..flow_solver.utils.boundary import cell_boundary as bdry_c
from ..flow_solver.physics import hydrostatics

from ..utils.data_structures import DiagnosticState


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

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = opts.BdryType.PERIODIC
        self.bdry_type[1] = opts.BdryType.WALL
        self.bdry_type[2] = opts.BdryType.PERIODIC

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

        self.diag_state = DiagnosticState(
            test_name="test_internal_long_wave",
            file_name="target_internal_long_wave",
            Nx=self.inx - 1,
            Ny=self.iny - 1,
            steps=[self.stepmax - 1],
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


def sol_init(Sol, npf, elem, node, th, ud, seeds=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed
    delth = 0.01 / ud.T_ref
    xc = -0.0 * ud.scale_factor * 50.0e3 / ud.h_ref
    xc = 0.0
    a = ud.scale_factor * 5.0e3 / ud.h_ref

    hydrostatics.analytical_state(npf, elem, node, th, ud)

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

    bdry_c.set_explicit_boundary_data(Sol, elem, ud, th, npf)

    return Sol


def T_from_p_rho(p, rho):
    return np.divide(p, rho)
