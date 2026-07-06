"""Williamson et al. (1992) test case 6 on the spherical channel.

Rossby-Haurwitz wave, wavenumber R = 4, as shallow water on the lat-lon
spherical channel (walls at +-80 degrees; the RH meridional velocity
scales as cos^3(phi) sin(phi), ~0.5% of the peak zonal speed at the
walls, so the channel adaptation perturbs the global solution weakly).
Same solver configuration as ``test_sphere_swe_tc2`` (thin frozen-radius
shell, tangent-plane constraint, traditional Coriolis f = 2 Omega
sin(phi)); see that module for the embedding/pseudovector conventions.

Nondimensionalization: h_ref = h0 = 8000 m, t_ref = sqrt(h0/g), so
g_swe = 1 and depth is O(1). The wave translates in longitude at the
analytic angular speed nu = [R(3+R)omega - 2 Omega]/[(1+R)(2+R)] — the
phase-speed gate of the validation scripts; the registered regression
run is SHORT (stepmax steps).
"""

import numpy as np

from ..utils import options as opts
from ..flow_solver.discretisation import spherical
from ..flow_solver.utils.boundary import node_boundary as bdry_n
from ..flow_solver.physics import hydrostatics

from .case_setup import build_bdry, do_initial_projection, make_diag_state

# --- Williamson TC6 constants (dimensional) -------------------------------
_G = 9.80616  # [m s^-2]
_A_EARTH = 6.37122e6  # [m]
_OMEGA_EARTH = 7.292e-5  # [s^-1]
_H0 = 8000.0  # [m]
_K = 7.848e-6  # [s^-1] (= omega, the wave angular parameters)
_R_WAVE = 4

# --- nondimensionalization -------------------------------------------------
T_REF = np.sqrt(_H0 / _G)
U_REF = _H0 / T_REF

A_ND = _A_EARTH / _H0
OMEGA_ND = _OMEGA_EARTH * T_REF
K_ND = _K * T_REF
PHI_MAX = np.deg2rad(80.0)
SHELL_HALF_DEPTH = 5.0

#: analytic zonal phase speed of the RH-4 pattern [rad / nondim time]
NU_PHASE = (_R_WAVE * (3 + _R_WAVE) * K_ND - 2.0 * OMEGA_ND) / (
    (1 + _R_WAVE) * (2 + _R_WAVE)
)


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

        self.g_swe = 1.0

        self.planet_radius = _A_EARTH
        self.curvilinear_map = spherical.SphericalShellMap(A_ND, frozen_radius=True)
        self.constrain_to_surface = True
        self.coriolis_field = self.curvilinear_map.traditional_coriolis(2.0 * OMEGA_ND)

        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.WALL
        )

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.45
        self.dtfixed = 2.0
        self.dtfixed0 = 2.0

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

        self.output_base_name = "_sphere_swe_tc6"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.diag_state = make_diag_state(
            "test_sphere_swe_tc6",
            "target_sphere_swe_tc6",
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


def rh4_fields(lam, phi):
    """Analytic RH-4 depth and (east, north) velocities, nondim,
    broadcastable. Standard Williamson TC6 formulas (their eq. 141-146)."""
    R = _R_WAVE
    om, K, Om, a = K_ND, K_ND, OMEGA_ND, A_ND
    cp = np.cos(phi)
    sp = np.sin(phi)

    u_east = a * om * cp + a * K * cp ** (R - 1) * (R * sp**2 - cp**2) * np.cos(R * lam)
    v_north = -a * K * R * cp ** (R - 1) * sp * np.sin(R * lam)

    A_t = om * (2.0 * Om + om) / 2.0 * cp**2 + K**2 / 4.0 * cp ** (2 * R) * (
        (R + 1) * cp**2 + (2 * R**2 - R - 2) - 2.0 * R**2 * cp ** (-2)
    )
    B_t = (
        2.0
        * (Om + om)
        * K
        / ((R + 1) * (R + 2))
        * cp**R
        * ((R**2 + 2 * R + 2) - (R + 1) ** 2 * cp**2)
    )
    C_t = K**2 / 4.0 * cp ** (2 * R) * ((R + 1) * cp**2 - (R + 2))

    # gh = gh0 + a^2 (A + B cos(R lam) + C cos(2R lam)); gh0 = 1 nondim
    gh = (
        1.0
        + a**2 * A_t
        + a**2 * B_t * np.cos(R * lam)
        + a**2 * C_t * np.cos(2 * R * lam)
    )
    return gh, u_east, v_north


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    g = ud.g_swe

    hydrostatics.integrated_state(npf, elem, node, th, ud)

    def coord(grid_obj, axis):
        shape = [1, 1, 1]
        shape[axis] = -1
        arr = (grid_obj.x, grid_obj.y, grid_obj.z)[axis]
        return arr.reshape(shape)

    lam, phi = coord(elem, 0), coord(elem, 2)
    gh, u_east, v_north = rh4_fields(lam, phi)
    h = gh / g

    # Cartesian momenta: east = e_lambda = (-sin lam, cos lam, 0),
    # north = e_phi = (-sin phi cos lam, -sin phi sin lam, -cos phi)
    # (velocities are true vectors; only pseudovectors flip in the
    # mirrored embedding)
    sl, cl = np.sin(lam), np.cos(lam)
    sp, cp = np.sin(phi), np.cos(phi)
    u_cart = (
        u_east * (-sl) + v_north * (-sp * cl),
        u_east * cl + v_north * (-sp * sl),
        v_north * (-cp) + 0.0 * lam,
    )

    shp = Sol.rho.shape
    Sol.rho[...] = np.broadcast_to(h, shp)
    Sol.rhou[...] = np.broadcast_to(h * u_cart[0], shp)
    Sol.rhov[...] = np.broadcast_to(h * u_cart[1], shp)
    Sol.rhow[...] = np.broadcast_to(h * u_cart[2], shp)
    Sol.rhoX[...] = 0.0

    p = g / 2.0 * Sol.rho**2
    Sol.rhoY[...] = p**th.gamminv

    lam_n, phi_n = coord(node, 0), coord(node, 2)
    gh_n, _, _ = rh4_fields(lam_n, phi_n)
    p_n = g / 2.0 * (gh_n / g) ** 2
    npf.p2_nodes[...] = np.broadcast_to(p_n**th.gamminv, npf.p2_nodes.shape)
    bdry_n.set_ghost_nodes(npf.p2_nodes, node, ud)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    do_initial_projection(Sol, npf, elem, node, th, ud, u0=0.0, v0=0.0)

    return Sol
