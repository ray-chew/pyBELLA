"""z-vertical plumbing smoke case (NOT a golden-master regression case).

A minimal 3D quasi-2D setup with the vertical on axis 2
(``gravity_direction = 2``, met convention): isothermal hydrostatic
background along z, a small thermal bump, a few steps. Its only job is to
push a non-default vertical axis through the PRODUCTION entry path —
``-ic`` registry, ``prepare.initialise`` (incl. ``axes.validate``),
hydrostatics along z, gravity ghost cells / wall zeroing on axis 2, and
the full time loop — which the in-process permutation oracle deliberately
bypasses. No target, no CompareSol: the pytest only asserts a clean run.

Physics agnosticity itself is proven by ``test_permutation_oracle.py``;
a permanent z-vertical golden-master case is deferred to the
terrain-following work (decision 2026-06-10).
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

        # vertical on axis 2 (z); x horizontal, y degenerate
        self.gravity_direction = 2
        self.gravity_strength = np.zeros(3)
        self.gravity_strength[2] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.i_gravity = np.zeros(3)
        self.i_gravity[2] = 1

        self.xmin = 0.0
        self.xmax = 60000.0 / self.h_ref
        self.ymin = 0.0
        self.ymax = 1.0
        self.zmin = 0.0
        self.zmax = 10000.0 / self.h_ref

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = opts.BdryType.PERIODIC
        self.bdry_type[1] = opts.BdryType.PERIODIC
        self.bdry_type[2] = opts.BdryType.WALL

        self.is_compressible = 1
        self.is_nonhydrostatic = 1

        self.CFL = 0.9
        self.dtfixed = 10.0 / self.t_ref  # 10 s
        self.dtfixed0 = 10.0 / self.t_ref

        self.inx = 48 + 1
        self.iny = 1 + 1
        self.inz = 10 + 1

        self.tout = [1.0]
        self.stepmax = 3

        self.stratification = self.stratification_function
        self.output_timesteps = False

        self.diag = False
        self.output_base_name = "_smoke_zvert"
        self.output_type = "test"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.inz - 1)
        self.autogen_fn = False

    def stratification_function(self, y):
        # isothermal: theta ~ exp(N^2 z / g) with N^2 = (gamma-1) g^2 / (gamma R T)
        Nsq = (
            ((self.gamm - 1.0) / self.gamm)
            * self.grav
            * self.grav
            / (self.R_gas * self.T_ref)
        ) * self.t_ref**2
        g = self.gravity_strength[2] / self.Msq
        return np.exp(Nsq * y / g)


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    # hydrostatic background along the configured vertical (z)
    hydrostatics.analytical_state(npf, elem, node, th, ud)

    v = axes.vertical_axis(ud)
    S0c = npf.HydroState.get_S0c(elem)
    rhoY0 = axes.expand_profile(npf.HydroState.rhoY0, elem.ndim, v, elem.sc)

    # small warm bump in the x-z plane
    x = elem.x.reshape(-1, 1, 1)
    z = elem.z.reshape(1, 1, -1)
    xc = 0.5 * (ud.xmax + ud.xmin)
    zc = 0.4 * ud.zmax
    delth = (0.1 / ud.T_ref) * np.exp(
        -((x - xc) ** 2 + (z - zc) ** 2) / (2000.0 / ud.h_ref) ** 2
    )

    Y = 1.0 / S0c + delth
    Sol.rhoY[...] = rhoY0
    Sol.rho[...] = rhoY0 / Y
    Sol.rhou[...] = 0.0
    Sol.rhov[...] = 0.0
    Sol.rhow[...] = 0.0
    Sol.rhoX[...] = Sol.rho * (Sol.rho / Sol.rhoY - S0c)

    npf.p2_nodes[...] = axes.expand_profile(npf.HydroState_n.p20, node.ndim, v, node.sc)

    ud.nonhydrostasy = 1.0
    ud.compressibility = 1.0

    return Sol
