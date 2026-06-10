"""Quasi-2D mountain-wave plumbing smoke case (NOT a golden-master case).

The Agnesi configuration *shape*: a vertical x-y slice run as quasi-2D 3D
(``inz = 2``, z degenerate periodic) so the implicit solve goes through the
full-tensor 27-point operator (``lap3D``) — the path terrain metric terms
attach to. Isothermal hydrostatic background, uniform horizontal wind,
periodic x, walls in y, a few steps. No target, no CompareSol: the pytest
only asserts a clean run.

Phase 0 of the terrain work runs it flat (no ``orography``); the witch-of-
Agnesi hill is switched on once the metric-aware operators land (Phase 6),
making this the first end-to-end terrain run. The golden-master Agnesi
regression case (``test_agnesi_hydrostatic``) is separate.
"""

import numpy as np

from ..utils import options as opts
from ..utils import axes
from ..flow_solver.physics import hydrostatics


class UserData(object):
    grav = 9.80665  # [m/s^2]
    omega = 0.0

    def __init__(self):
        self.grav = self.grav
        self.omega = self.omega

        self.R_gas = 287.05
        self.gamm = 1.4
        self.T_ref = 250.0
        self.p_ref = 1e5
        self.h_ref = self.R_gas * self.T_ref / self.grav  # [m]
        self.u_ref = 10.0
        self.t_ref = self.h_ref / self.u_ref

        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        # vertical on axis 1 (y); x horizontal, z degenerate
        self.gravity_direction = 1
        self.gravity_strength = np.zeros(3)
        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.i_gravity = np.zeros(3)
        self.i_gravity[1] = 1

        self.xmin = -30000.0 / self.h_ref
        self.xmax = 30000.0 / self.h_ref
        self.ymin = 0.0
        self.ymax = 10000.0 / self.h_ref
        self.zmin = 0.0
        self.zmax = 1.0

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = opts.BdryType.PERIODIC
        self.bdry_type[1] = opts.BdryType.WALL
        self.bdry_type[2] = opts.BdryType.PERIODIC

        self.is_compressible = 1
        self.is_nonhydrostatic = 1

        # uniform background wind along x
        self.u_wind_speed = 10.0 / self.u_ref
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.CFL = 0.9
        self.dtfixed = 10.0 / self.t_ref  # 10 s
        self.dtfixed0 = 10.0 / self.t_ref

        self.inx = 48 + 1
        self.iny = 16 + 1
        self.inz = 1 + 1

        self.tout = [1.0]
        self.stepmax = 3

        self.stratification = self.stratification_function
        self.output_timesteps = False

        self.diag = False
        self.output_base_name = "_smoke_agnesi"
        self.output_type = "test"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)
        self.autogen_fn = False

    def stratification_function(self, y):
        # isothermal: theta ~ exp(N^2 z / g) with N^2 = (gamma-1) g^2 / (gamma R T)
        Nsq = (
            ((self.gamm - 1.0) / self.gamm)
            * self.grav
            * self.grav
            / (self.R_gas * self.T_ref)
        ) * self.t_ref**2
        g = self.gravity_strength[1] / self.Msq
        return np.exp(Nsq * y / g)


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    # hydrostatic background along the configured vertical (y)
    hydrostatics.analytical_state(npf, elem, node, th, ud)

    v = axes.vertical_axis(ud)
    S0c = npf.HydroState.get_S0c(elem)
    rhoY0 = axes.expand_profile(npf.HydroState.rhoY0, elem.ndim, v, elem.sc)

    Y = 1.0 / S0c
    Sol.rhoY[...] = rhoY0
    Sol.rho[...] = rhoY0 / Y
    Sol.rhou[...] = Sol.rho * ud.u_wind_speed
    Sol.rhov[...] = 0.0
    Sol.rhow[...] = 0.0
    Sol.rhoX[...] = Sol.rho * (Sol.rho / Sol.rhoY - S0c)

    npf.p2_nodes[...] = axes.expand_profile(npf.HydroState_n.p20, node.ndim, v, node.sc)

    ud.nonhydrostasy = 1.0
    ud.compressibility = 1.0

    return Sol
