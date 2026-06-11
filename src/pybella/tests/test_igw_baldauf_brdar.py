"""Baldauf & Brdar (2013) linear internal gravity wave test.

Ported from the legacy ``igw_baldauf_brdar`` initial condition
(recoverable from the git tag ``archive/full_coriolis``). The case is a
small thermal perturbation on an isothermal (constant Brunt-Väisälä
frequency) hydrostatic background in a periodic channel, for which
Baldauf & Brdar (QJRMS, 2013) derive an analytic linear solution.

Physics deltas relative to ``test_internal_long_wave``:

- isothermal background: ``Nsq_ref = ((gamma-1)/gamma) * g^2 / (R_gas * T_ref)``
  with ``T_ref = 250 K`` (instead of the prescribed ``Nsq_ref = 1e-4``),
- reference scales ``u_ref = 10 m/s``, ``h_ref = R_gas * T_ref / grav``,
  ``t_ref = h_ref / u_ref``,
- domain: ``scale_factor * 300 km`` wide, ``10 km`` deep (long-wave
  variant for ``scale_factor = 20``),
- Gaussian envelope perturbation
  ``delT * molly(x) * sin(pi*y/H) * exp(-(x-xc)^2/a^2)`` instead of the
  Lorentzian ``1/(1+(x-xc)^2/a^2)``,
- compressible (``is_compressible = 1``),
- traditional-Coriolis component ``coriolis_strength[1] = omega * t_ref``
  with ``omega = 1.03126e-4 s^-1``.
"""

import numpy as np

from ..utils import options as opts

from ..flow_solver.utils import fields
from ..flow_solver.physics import hydrostatics

from ..utils.data_structures import DiagnosticState


class UserData(object):
    # planetary -> 160.0;  long-wave -> 20.0;  standard -> 1.0;
    scale_factor = 20.0

    def __init__(self):
        self.scale_factor = self.scale_factor

        self.grav = 9.80665  # [m/s^2]
        self.omega = 1.0 * 0.000103126  # [s^{-1}]
        self.R_gas = 287.05  # [J kg^{-1} K^{-1}]
        self.gamm = 1.4

        self.T_ref = 250.00  # [K]
        self.u_ref = 10.0  # [m/s]
        self.p_ref = 1e5  # [Pa]
        self.h_ref = self.R_gas * self.T_ref / self.grav  # [m]
        self.t_ref = self.h_ref / self.u_ref  # [s]

        self.Nsq_ref = (
            ((self.gamm - 1.0) / self.gamm)
            * self.grav
            * self.grav
            / (self.R_gas * self.T_ref)
        )  # [s^{-2}]

        self.Msq = (
            self.u_ref * self.u_ref / (self.R_gas * self.T_ref)
        )  # Mach number squared

        self.is_nonhydrostatic = 1
        self.is_compressible = 1

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

        # NB: defined after `omega` and `t_ref` so that the dependency
        # manager in `UserDataInit` does not overwrite this explicit choice.
        self.coriolis_strength = np.zeros((3))
        self.coriolis_strength[1] = self.omega * self.t_ref

        self.xmin = -0.5 * self.scale_factor * 300000.0 / self.h_ref
        self.xmax = 0.5 * self.scale_factor * 300000.0 / self.h_ref
        self.ymin = 0.0
        self.ymax = 10000.0 / self.h_ref
        self.zmin = -1.0
        self.zmax = 1.0

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = opts.BdryType.PERIODIC
        self.bdry_type[1] = opts.BdryType.WALL
        self.bdry_type[2] = opts.BdryType.WALL

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9

        self.dtfixed0 = 0.5 / self.t_ref * 50.0 * self.scale_factor
        self.dtfixed = 0.5 / self.t_ref * 50.0 * self.scale_factor

        self.inx = 301 + 1
        self.iny = 20 + 1
        self.inz = 1

        self.tout = [8 * 180.0 * self.scale_factor / self.t_ref]  # 8 hrs

        self.tol = 1.0e-8
        self.stepmax = 31
        self.max_iterations = 6000

        self.autogen_fn = False

        self.output_timesteps = True

        self.stratification = self.stratification_function
        self.molly = self.molly_function
        self.rhoe = self.rhoe_method

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_igw_baldauf_brdar"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.diag_state = DiagnosticState(
            test_name="test_igw_baldauf_brdar",
            file_name="target_igw_baldauf_brdar",
            Nx=self.inx - 1,
            Ny=self.iny - 1,
            steps=[self.stepmax - 1],
            # tolerance audit 2026-06-09 measured ~1e-10 run-to-run scatter
            # SAME-MACHINE; the first CI run (2026-06-11, GitHub runner)
            # showed cross-PLATFORM scatter of 2.3e-6 on rhou — different
            # CPU/BLAS/numba reorder the bicgstab reductions. Gate at the
            # 1e-5 default; physics is guarded by the B&B analytic oracle.
            tolerances={
                "rho": 1e-5,
                "rhou": 1e-5,
                "rhov": 1e-5,
                "rhow": 1e-5,
                "rhoY": 1e-5,
                "rhoX": 1e-5,
                "p2_nodes": 1e-5,
            },
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
    delT = 0.01 / ud.T_ref

    xc = 0.0
    a = ud.scale_factor * 5.0e3 / ud.h_ref
    H = ud.ymax - ud.ymin

    hydrostatics.analytical_state(npf, elem, node, th, ud)

    x = elem.x.reshape(-1, 1)
    y = elem.y.reshape(1, -1)

    Tb = delT * ud.molly(x) * np.sin(np.pi * y / H) * np.exp(-((x - xc) ** 2) / (a**2))

    xn = node.x[:-1].reshape(-1, 1)
    yn = node.y[:-1].reshape(1, -1)

    Tbn = (
        delT
        * ud.molly(xn)
        * np.sin(np.pi * yn / H)
        * np.exp(-((xn - xc) ** 2) / (a**2))
    )

    HySt = fields.States(node.sc)
    HyStn = fields.States(node.sc)

    Y = ud.stratification(y) + Tb
    Yn = ud.stratification(yn) + Tbn

    hydrostatics.column(HySt, HyStn, Y, Yn, elem, node, th, ud)

    xc_idx = slice(0, -1)
    yc_idx = slice(0, -1)
    c_idx = (xc_idx, yc_idx)

    if ud.is_compressible:
        p = HySt.p0[c_idx]
        rhoY = HySt.rhoY0[c_idx]
    else:
        p = npf.HydroState.p0
        rhoY = npf.HydroState.rhoY0

    u, v, w = u0, v0, w0

    rho = rhoY / Y
    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhow[...] = rho * w
    Sol.rhoY[...] = rhoY

    npf.p2_cells[...] = HySt.p20[c_idx]

    Sol.rhoX[...] = Sol.rho * (Sol.rho / Sol.rhoY - npf.HydroState.S0.reshape(1, -1))

    npf.p2_nodes[:, :] = HyStn.p20

    hydrostatics.initial_pressure(Sol, npf, elem, node, ud, th)

    ud.nonhydrostasy = 1.0 if ud.is_nonhydrostatic == 1 else 0.0
    ud.compressibility = 1.0 if ud.is_compressible == 1 else 0.0

    if "imbal" in ud.aux:
        npf.p2_nodes[...] = 0.0

    return Sol


def T_from_p_rho(p, rho):
    return np.divide(p, rho)
