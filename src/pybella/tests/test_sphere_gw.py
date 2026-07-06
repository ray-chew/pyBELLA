"""DCMIP-31-style nonhydrostatic gravity wave on the small-planet shell.

The 3D compressible spherical shell (``SphericalShellMap``, radial gravity)
carries a small potential-temperature perturbation on the isothermal
(constant Brunt-Vaisala) hydrostatic background; the imbalance launches a
gravity-wave train that propagates zonally around the reduced-radius planet
(DCMIP 2012 test 3-1 geometry: a = a_earth / X, X = 125, no rotation).

The perturbation is deliberately LATITUDE-INDEPENDENT and periodic in
longitude,

    theta'(lambda, r) = delTheta * exp(-(lambda_wrapped / Lam)^2)
                                 * sin(pi * (r - a) / depth),

so that in the large-radius limit (a -> inf, equatorial band) the near-
equator lambda-height slice reduces to the planar Baldauf & Brdar (2013)
isothermal channel with x = a*lambda -- the physics oracle for the
``test_igw_baldauf_brdar`` machinery (``baldauf_brdar_analytic``). The
thermodynamic constants are shared with that case (T_ref = 250 K,
u_ref = 10 m/s, h_ref = R T_ref / g) so the two backgrounds coincide.

Nondimensionalization (isothermal scale height): h_ref = R_gas*T_ref/grav,
t_ref = h_ref/u_ref, gravity_strength[1] = g*h_ref/(R*T_ref) = 1. The shell
is the true r-dependent metric (``frozen_radius=False``); gravity stays
radial (``up_direction`` = e_r). Balanced background => p2 = 0 at init and
zero initial velocity; the wave grows from the theta perturbation alone.

Overridable knobs (mutated by the in-process gates in
``test_scripts/test_sphere_gw.py``): ``X`` (radius reduction factor),
``pert_amplitude`` (delTheta / T_ref), the grid counts and ``stepmax``.
"""

import numpy as np

from ..utils import options as opts
from ..flow_solver.discretisation import spherical
from ..flow_solver.physics import hydrostatics

from .case_setup import build_bdry, make_diag_state

_A_EARTH = 6.37122e6  # [m]


class UserData(object):
    grav = 9.80665  # [m/s^2]
    omega = 0.0  # nonrotating (DCMIP 3-1)

    #: planet-radius reduction factor a = a_earth / X (DCMIP 3-1: 125)
    X = 125.0
    #: potential-temperature perturbation amplitude delTheta / T_ref
    pert_amplitude = 0.01 / 250.0
    #: zonal angular half-width of the perturbation [rad]
    pert_halfwidth = 0.4
    #: shell depth (radial extent) [m]
    depth_m = 10000.0
    #: latitude band half-extent [rad] (equatorial channel; walls in phi)
    phi_band = 0.3

    def __init__(self):
        self.grav = self.grav
        self.omega = self.omega
        self.X = self.X
        self.pert_amplitude = self.pert_amplitude
        self.pert_halfwidth = self.pert_halfwidth
        self.depth_m = self.depth_m
        self.phi_band = self.phi_band

        self.R_gas = 287.05  # [J kg^-1 K^-1]
        self.gamm = 1.4
        self.T_ref = 250.0  # [K]
        self.u_ref = 10.0  # [m/s]
        self.p_ref = 1e5  # [Pa]
        self.h_ref = self.R_gas * self.T_ref / self.grav  # [m]
        self.t_ref = self.h_ref / self.u_ref  # [s]

        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.is_compressible = 1
        self.is_nonhydrostatic = 1

        # radial gravity: gravity_direction = 1 (axis y = r), unit magnitude
        self.gravity_direction = 1
        self.gravity_strength = np.zeros(3)
        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.i_gravity = np.zeros(3)
        self.i_gravity[1] = 1

        self.planet_radius = _A_EARTH / self.X
        a_nd = self.planet_radius / self.h_ref
        depth = self.depth_m / self.h_ref

        # axes: x = lambda, y = r, z = phi
        self.xmin, self.xmax = -np.pi, np.pi
        self.ymin, self.ymax = a_nd, a_nd + depth
        self.zmin, self.zmax = -self.phi_band, self.phi_band

        self.curvilinear_map = spherical.SphericalShellMap(a_nd, frozen_radius=False)

        self.u_wind_speed = 0.0
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.WALL
        )

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.45
        # horizontal acoustic CFL on the small planet sets dt; calibrated so
        # a 64x8x8 run is comfortably stable, scaled with the grid by the gates
        self.dtfixed = 5.0 / self.t_ref
        self.dtfixed0 = 5.0 / self.t_ref

        self.inx = 64 + 1
        self.iny = 8 + 1
        self.inz = 8 + 1

        self.initial_projection = False

        self.tout = [1.0e9]
        self.stepmax = 60

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_function
        self.output_timesteps = True

        self.diag = False
        self.diag_updt_targets = False

        self.output_base_name = "_sphere_gw"
        self.output_type = "test"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.inz - 1)
        self.autogen_fn = False

        self.diag_state = make_diag_state(
            "test_sphere_gw",
            "target_sphere_gw",
            self.inx,
            self.iny,
            self.stepmax,
        )

    def stratification_function(self, z):
        # isothermal: theta ~ exp(N^2 z / g), constant Brunt-Vaisala
        Nsq = (
            ((self.gamm - 1.0) / self.gamm)
            * self.grav
            * self.grav
            / (self.R_gas * self.T_ref)
        ) * self.t_ref**2
        g = self.gravity_strength[1] / self.Msq
        return np.exp(Nsq * z / g)

    def rhoe_function(self, rho, u, v, w, p, ud, th):
        Msq = ud.compressibility * ud.Msq
        gm1inv = th.gm1inv
        return p * gm1inv + 0.5 * Msq * rho * (u**2 + v**2 + w**2)


def _theta_perturbation(elem, ud):
    """Latitude-independent, longitude-periodic theta' field (cell grid)."""
    a_nd = ud.curvilinear_map.radius
    depth = ud.ymax - ud.ymin
    height = elem.metric.height  # r - a, cell field (elem.sc)

    lam = elem.x.reshape(-1, 1, 1)
    lam_w = np.mod(lam + np.pi, 2.0 * np.pi) - np.pi

    envelope = np.exp(-((lam_w / ud.pert_halfwidth) ** 2))
    vertical = np.sin(np.pi * height / depth)
    return ud.pert_amplitude * envelope * vertical


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    hydrostatics.analytical_state(npf, elem, node, th, ud)

    S0c = npf.HydroState.get_S0c(elem)  # background inverse-theta (field)
    rhoY0 = npf.HydroState.rhoY0  # hydrostatic pressure variable (field)
    Y0 = npf.HydroState.Y0  # background potential temperature (field)

    theta_pert = _theta_perturbation(elem, ud)
    Y = Y0 + theta_pert

    Sol.rhoY[...] = rhoY0
    Sol.rho[...] = rhoY0 / Y  # perturb theta at constant pressure
    Sol.rhou[...] = 0.0
    Sol.rhov[...] = 0.0
    Sol.rhow[...] = 0.0
    Sol.rhoX[...] = Sol.rho * (Sol.rho / Sol.rhoY - S0c)

    npf.p2_nodes[...] = 0.0

    ud.nonhydrostasy = 1.0 if ud.is_nonhydrostatic else 0.0
    ud.compressibility = 1.0 if ud.is_compressible else 0.0

    return Sol
