"""Williamson et al. (1992) test case 2 on the spherical channel.

Steady-state geostrophically balanced zonal flow, run as shallow water on
the lat-lon spherical channel (all longitudes, walls at +-80 degrees
latitude) through the gamma = 2 SWE equivalence of ``test_swe_vortex``::

    rho  <->  fluid depth h        p    <->  g h^2 / 2
    rhoY <->  p^(1/gamma) = sqrt(g/2) h

The sphere enters as metric data (``SphericalShellMap``, thin frozen-
radius shell, one interior radial cell) with Cartesian momenta constrained
to the local tangent plane (``ud.constrain_to_surface``) and the
traditional-approximation Coriolis field f(phi) e_r, f = 2 Omega sin(phi).

Nondimensionalization (values pre-scaled; the ud refs stay 1 so Msq = 1):
h_ref = H0 = gh0/g = 2998.11 m, t_ref = sqrt(H0/g) = 17.485 s,
u_ref = sqrt(g H0) = 171.47 m/s. Then g_swe = 1, gh0 = 1,
a = 2124.9, Omega = 1.27506e-3, u0 = 2 pi a / 12 days = 0.22517.

The analytic steady state is

    gh(phi) = gh0 - (a Omega u0 + u0^2/2) sin^2(phi),
    u = u0 cos(phi) e_lambda   (Cartesian components; e_lambda = (-sin
    lambda, cos lambda, 0) in the pole-along-minus-x2 embedding).

The registered regression run is SHORT (stepmax steps); the 12-day
l2(h) <= 1e-3 validation gate runs from ``test_scripts`` /
``run_scripts`` where wall-clock allows.
"""

import numpy as np

from ..utils import options as opts
from ..flow_solver.discretisation import spherical
from ..flow_solver.utils.boundary import node_boundary as bdry_n
from ..flow_solver.physics import hydrostatics

from .case_setup import build_bdry, do_initial_projection, make_diag_state

# --- Williamson TC2 constants (dimensional) -------------------------------
_G = 9.80616  # [m s^-2]
_A_EARTH = 6.37122e6  # [m]
_OMEGA_EARTH = 7.292e-5  # [s^-1]
_GH0 = 2.94e4  # [m^2 s^-2]

# --- nondimensionalization -------------------------------------------------
H0 = _GH0 / _G  # h_ref [m]
T_REF = np.sqrt(H0 / _G)  # [s]
U_REF = H0 / T_REF  # = sqrt(g H0) [m/s]

A_ND = _A_EARTH / H0
OMEGA_ND = _OMEGA_EARTH * T_REF
U0_ND = (2.0 * np.pi * _A_EARTH / (12.0 * 86400.0)) / U_REF
PHI_MAX = np.deg2rad(80.0)
SHELL_HALF_DEPTH = 10.0  # nondim radial half-thickness of the shell


class UserData(object):
    grav = 0.0
    omega = 0.0

    R_gas = 1.0
    gamm = 2.0

    h_ref = 1.0
    t_ref = 1.0
    T_ref = 1.0
    p_ref = 1.0

    def __init__(self):
        self.grav = self.grav
        self.omega = self.omega
        self.R_gas = self.R_gas
        self.gamm = self.gamm
        self.h_ref = self.h_ref
        self.t_ref = self.t_ref
        self.T_ref = self.T_ref
        self.p_ref = self.p_ref

        # lambda, r, phi (axes x, y, z)
        self.xmin = -np.pi
        self.xmax = np.pi
        self.ymin = A_ND - SHELL_HALF_DEPTH
        self.ymax = A_ND + SHELL_HALF_DEPTH
        self.zmin = -PHI_MAX
        self.zmax = PHI_MAX

        self.u_wind_speed = 0.0
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        # SWE gravitational acceleration (nondim; pressure law only)
        self.g_swe = 1.0

        self.planet_radius = _A_EARTH
        self.curvilinear_map = spherical.SphericalShellMap(A_ND, frozen_radius=True)
        self.constrain_to_surface = True
        # traditional approximation, kernel convention w = f = 2 Omega
        self.coriolis_field = self.curvilinear_map.traditional_coriolis(2.0 * OMEGA_ND)

        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.WALL
        )

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.45
        self.dtfixed = 5.0
        self.dtfixed0 = 5.0

        self.inx = 128 + 1
        self.iny = 2  # one interior radial cell (thin shell)
        self.inz = 64 + 1

        self.initial_projection = True

        self.tout = [1.0e9]
        self.stepmax = 31

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_function
        self.output_timesteps = True

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_sphere_swe_tc2"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.diag_state = make_diag_state(
            "test_sphere_swe_tc2",
            "target_sphere_swe_tc2",
            self.inx,
            self.iny,
            self.stepmax,
            plot_compare=True,
        )

        self.autogen_fn = False

    def stratification_function(self, y):
        if type(y) == float:
            return 1.0
        else:
            return np.ones((y.shape))

    def rhoe_function(self, rho, u, v, w, p, ud, th):
        Msq = ud.compressibility * ud.Msq
        gm1inv = th.gm1inv

        return p * gm1inv + 0.5 * Msq * rho * (u**2 + v**2 + w**2)


def depth_and_velocity(lam, phi):
    """Analytic TC2 fields (nondim): depth h(phi) and the Cartesian
    momentum direction of u0 cos(phi) e_lambda. Broadcastable."""
    c_h = A_ND * OMEGA_ND * U0_ND + 0.5 * U0_ND**2
    h = 1.0 - c_h * np.sin(phi) ** 2
    u_cart = (
        U0_ND * np.cos(phi) * (-np.sin(lam)),
        U0_ND * np.cos(phi) * np.cos(lam) + 0.0 * phi,
        0.0 * (lam + phi),
    )
    return h, u_cart


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    g = ud.g_swe

    hydrostatics.integrated_state(npf, elem, node, th, ud)

    def coord(grid_obj, axis):
        shape = [1, 1, 1]
        shape[axis] = -1
        arr = (grid_obj.x, grid_obj.y, grid_obj.z)[axis]
        return arr.reshape(shape)

    # cells: fill EVERYWHERE analytically (lambda ghosts wrap through the
    # trig functions; phi/r ghosts are overwritten by set_ghost_cells)
    lam_c, phi_c = coord(elem, 0), coord(elem, 2)
    h, u_cart = depth_and_velocity(lam_c, phi_c)

    shp = Sol.rho.shape
    Sol.rho[...] = np.broadcast_to(h, shp)
    Sol.rhou[...] = np.broadcast_to(h * u_cart[0], shp)
    Sol.rhov[...] = np.broadcast_to(h * u_cart[1], shp)
    Sol.rhow[...] = np.broadcast_to(h * u_cart[2], shp)
    Sol.rhoX[...] = 0.0

    p = g / 2.0 * Sol.rho**2
    Sol.rhoY[...] = p**th.gamminv

    # nodes: same analytic depth on the node grid
    lam_n, phi_n = coord(node, 0), coord(node, 2)
    h_n, _ = depth_and_velocity(lam_n, phi_n)
    p_n = g / 2.0 * h_n**2
    npf.p2_nodes[...] = np.broadcast_to(p_n**th.gamminv, npf.p2_nodes.shape)
    bdry_n.set_ghost_nodes(npf.p2_nodes, node, ud)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    do_initial_projection(Sol, npf, elem, node, th, ud, u0=0.0, v0=0.0)

    return Sol
