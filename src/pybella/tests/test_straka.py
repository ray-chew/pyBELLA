"""Straka density current (Straka et al. 1993).

Cold-bubble benchmark: a -15 K temperature anomaly dropped into a neutrally
stratified (theta = 300 K), hydrostatically balanced atmosphere collapses,
spreads along the ground and rolls up into Kelvin-Helmholtz rotors. The
benchmark prescribes a fixed physical viscosity/diffusivity K = 75 m^2/s on
velocity and potential temperature so the solution converges with
resolution; this case is the reason ``flow_solver/numerics/diffusion.py``
exists (``ud.diffusion`` flag).

Setup (dimensional): domain 51.2 km x 6.4 km, free-slip walls all around,
dx = dz = 200 m, dt = 4 s fixed, run to t = 900 s. Anomaly centred at
x = 0 km, z = 3 km with radii (4 km, 2 km):

    T' = -15/2 * (1 + cos(pi * r))  [K]   for  r <= 1,
    theta' = T' / pi_bar(z)               (perturbation at fixed pressure).

Non-dimensionalisation: h_ref = 10 km, t_ref = 1000 s (u_ref = 10 m/s),
T_ref = 300 K, p_ref = 1e5 Pa; K* = K t_ref / h_ref^2 = 7.5e-4.

This is the suite's only nonlinear, advection-dominated gravity+wall case.
It supersedes the ``rising_bubble_cold`` IC (recoverable from the git tag
``archive/full_coriolis``).
"""

import numpy as np

from ..flow_solver.physics import hydrostatics
from ..utils import options as opts
from .case_setup import build_bdry, make_diag_state


class UserData(object):
    def __init__(self):
        self.grav = 9.81  # [m/s^2]
        self.t_ref = 1000.0  # [s]
        self.T_ref = 300.0  # [K]
        self.h_ref = 10000.0  # [m]
        self.p_ref = 1e5  # [Pa]

        self.xmin = -2.56  # [-25.6 km]
        self.xmax = 2.56
        self.ymin = 0.0
        self.ymax = 0.64  # [6.4 km]

        # free-slip walls all around, faithful to Straka et al. (1993).
        # x-WALLs were broken until the axial-agnosticity boundary fixes
        # (the nodal-divergence wall zeroing was vertical-axis-only and the
        # wall-normal momentum mirror was hardcoded to rhov); this case now
        # exercises that path. With the domain at +-25.6 km the fronts
        # (~15.5 km at t=900s) stay well clear of the boundary either way.
        self.bdry_type = build_bdry(
            opts.BdryType.WALL, opts.BdryType.WALL, opts.BdryType.WALL
        )

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9  # cap only; dtfixed binds for |u| < 45 m/s
        self.dtfixed = 0.004  # 4 s
        self.dtfixed0 = 0.004

        self.inx = 256 + 1  # dx = 200 m
        self.iny = 32 + 1  # dy = 200 m
        self.inz = 1

        self.tout = [0.9]  # 900 s
        # exactly 225 steps (000..224); strip_target_file keeps step stepmax-1
        self.stepmax = 225

        self.is_compressible = 1

        # Straka's fixed physical viscosity/diffusivity K = 75 m^2/s
        self.diffusion = True
        self.diffusion_coeff = 75.0 * self.t_ref / self.h_ref**2  # 7.5e-4

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_straka"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.output_timesteps = True

        self.diag_state = make_diag_state(
            "test_straka",
            "target_straka",
            self.inx,
            self.iny,
            self.stepmax,
        )

        self.autogen_fn = False


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed

    delT = -15.0  # [K] temperature anomaly amplitude
    xc, yc = 0.0, 0.3  # centre: (0 km, 3 km)
    xr, yr = 0.4, 0.2  # radii: (4 km, 2 km)

    hydrostatics.integrated_state(npf, elem, node, th, ud)

    x = elem.x
    y = elem.y
    x, y = np.meshgrid(x, y)

    r = np.sqrt(((x - xc) / xr) ** 2 + ((y - yc) / yr) ** 2)

    # temperature anomaly at fixed pressure -> theta' = T'/pi_bar
    perturbation = (delT / ud.T_ref) * 0.5 * (np.cos(np.pi * r) + 1.0)
    perturbation[np.where(r > 1.0)] = 0.0

    rhoY = npf.HydroState.rhoY0[np.newaxis, :]
    pi_bar = rhoY**th.gm1

    rho = rhoY / (ud.stratification(y) + perturbation.T / pi_bar)

    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u0
    Sol.rhov[...] = rho * v0
    Sol.rhow[...] = rho * w0
    Sol.rhoY[...] = rhoY

    # hydrostatically balanced background: zero perturbation pressure
    npf.p2_nodes[...] = 0.0

    return Sol
