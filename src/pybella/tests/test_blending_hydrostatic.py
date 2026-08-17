"""Hydrostatic <-> nonhydrostatic blending regression case (thesis sec. 4.2).

Imbalanced inertia-gravity wave in the Skamarock & Klemp (1994) *hydrostatic*
configuration (wide aspect ratio), the demonstration archetype of thesis
sec. 6.2.3. Physics is the ``test_internal_long_wave`` long-wave IGW; the
blending setup mirrors ``test_blending_warm_bubble``:

* the initial pressure is left imbalanced w.r.t. the potential-temperature
  perturbation (``p2_nodes = 0``, triggered by ``"imbal"`` in ``ud.aux``);
* the schedule-driven blend (``continuous_blending = True``,
  ``no_of_hy_initial = 1``) runs the first step in the hydrostatic regime
  (``is_nonhydrostatic = 0``), then flips to nonhydrostatic -- the hydrostatic
  elliptic operator + explicit vertical-momentum switch carry the balance, with
  NO explicit conversion routine (the mechanism inherited from the predecessor
  RKLM_Python code);
* ``no_of_pi_initial = 0`` disables comp<->psinc blending so only the
  hydro<->nonhydro switch is exercised.

The hydrostatic background MUST use ``hydrostatics.integrated_state`` (the
stratification-consistent quadrature), not the isothermal ``analytical_state`` --
this case is constant-N. The isothermal background gave N^2 3.19x too large and a
huge spurious ``rhoX`` that broke the blend.

Like the warm-bubble case this is a *reproducibility gate*: all dynamic fields
are un-gated (tol 1.0) and only ``p2_nodes`` is gated (the blending must keep the
vertically-propagating acoustic mode suppressed), compared as a time increment.
"""

import numpy as np

from ..utils import options as opts

from ..flow_solver.utils import fields
from ..flow_solver.physics import hydrostatics

from .case_setup import build_bdry, make_diag_state


class UserData(object):
    # planetary -> 160.0;  long-wave / hydrostatic -> 20.0;  standard -> 1.0;
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
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.gravity_strength = np.zeros((3))
        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)

        gravity_mask = (self.gravity_strength > np.finfo(np.float64).eps) | (
            np.arange(3) == 1
        )
        self.i_gravity = gravity_mask.astype(int)
        if np.any(gravity_mask):
            self.gravity_direction = np.where(gravity_mask)[0][-1]

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
        self.dtfixed0 = 1.0
        self.dtfixed = 1.0

        self.inx = 301 + 1
        self.iny = 10 + 1
        self.inz = 1

        self.tout = [self.scale_factor * 1.0 * 3000.0 / self.t_ref]

        self.tol = 1.0e-12
        self.stepmax = 11
        self.max_iterations = 6000

        ##########################################
        # BLENDING : hydrostatic -> nonhydrostatic
        ##########################################
        self.is_compressible = 1
        self.is_nonhydrostatic = 1

        # Schedule-driven blend (inherited from RKLM_Python): the eos
        # schedule (physics/eos.py) sets is_nonhydrostatic = 0 for the first
        # ``no_of_hy_initial`` steps then flips to 1 -- the hydrostatic elliptic
        # operator + explicit vertical-momentum switch carry the balance, with
        # NO explicit hydro<->nonhydro conversion routine. Requires
        # continuous_blending = True; initial_blending stays False.
        self.continuous_blending = True
        self.no_of_pi_initial = 0  # no comp<->psinc blending
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 1  # one hydrostatic step, then blend to nonhydro
        self.no_of_hy_transition = 0

        self.initial_blending = False

        self.autogen_fn = False
        self.output_timesteps = True

        self.stratification = self.stratification_function
        self.molly = self.molly_function
        self.rhoe = self.rhoe_method

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_blending_hydrostatic"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        # "imbal": start from an imbalanced pressure field and route step 0
        # through the initial hydrostatic projection (see schemes/orchestration).
        self.aux = "imbal"
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.diag_state = make_diag_state(
            f"{self.output_type}_blending_hydrostatic",
            "target_blending_hydrostatic",
            self.inx,
            self.iny,
            self.stepmax,
            plot_compare=True,
            time_increment=True,
            # The only thing that matters here is that p2_nodes remains small,
            # i.e. the blend suppresses the vertically-propagating acoustic mode.
            tolerances={
                "rho": 1.0e-0,
                "rhou": 1.0e-0,
                "rhov": 1.0e-0,
                "rhow": 1.0e-0,
                "rhoY": 1.0e-0,
                "rhoX": 1.0e-0,
                "p2_nodes": 1.0e-4,
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


def sol_init(Sol, npf, elem, node, th, ud, seeds=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed
    delth = 0.01 / ud.T_ref
    xc = 0.0
    a = ud.scale_factor * 5.0e3 / ud.h_ref

    hydrostatics.integrated_state(npf, elem, node, th, ud)

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

    # imbalanced initial data: the pressure is NOT balanced w.r.t. the
    # potential-temperature perturbation (thesis sec. 6.2). The hydrostatic
    # blending step is what recovers a balanced pressure/vertical momentum.
    if "imbal" in ud.aux:
        npf.p2_nodes[...] = 0.0

    return Sol
