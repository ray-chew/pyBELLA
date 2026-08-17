"""Schär ridge (Schär et al. 2002) — SLEVE golden-master case, native 2D.

The classic two-scale mountain-wave benchmark: a Gaussian envelope
modulated by a small-scale cosine ridge,

    h(x) = h0 exp(-(x/a)^2) cos^2(pi x / lambda),

h0 = 250 m, a = 5 km, lambda = 4 km, in uniform wind U = 10 m/s with
constant buoyancy frequency N = 0.01 1/s. Spectrally the terrain has an
envelope peak at k ~ 0 (hydrostatic, radiates) and small-scale peaks at
k = +-2 pi / lambda — EVANESCENT here (k2 = 1.57e-3 > l = N/U = 1e-3
1/m), so above ~2 km the true solution carries no small-scale signal.
That is what discriminates vertical coordinates: under Gal-Chen the
small-scale terrain distorts the grid at every level and leaves spurious
small-scale w aloft; under SLEVE (this case's transform) the distortion
decays on s2 = 2.5 km and the upper levels stay clean.

The exact cos^2 identity supplies the SLEVE split in closed form:

    smooth   h1(x) = (h0/2) exp(-(x/a)^2)
    residual h2(x) = (h0/2) exp(-(x/a)^2) cos(2 pi x / lambda)

with analytic gradients for both (max slope ~0.2 — 25x the Agnesi case;
FD slopes are not trusted at that steepness).

Native 2D (lap2D terrain path); the regression run is a determinism gate.
Physics + the Gal-Chen-vs-SLEVE discriminator live in
``test_scripts/test_schaer_analytic.py`` against the linear FFT oracle
(``tests/schaer_linear_analytic.py``).
"""

import numpy as np

from ..utils import options as opts
from ..flow_solver.physics import hydrostatics
from ..flow_solver.discretisation import terrain

from .case_setup import apply_rayleigh_bdry, build_bdry, make_diag_state

_COMPARED_FIELDS = ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX", "p2_nodes")


class UserData(object):
    grav = 9.81  # [m/s^2]
    omega = 0.0

    NN = 0.01  # [1/s] buoyancy frequency
    U0 = 10.0  # [m/s] background wind
    hill_height = 250.0  # [m] h0
    hill_width = 5000.0  # [m] envelope half-width a
    ridge_wavelength = 4000.0  # [m] lambda

    def __init__(self):
        self.grav = self.grav
        self.omega = self.omega
        self.NN = self.NN
        self.U0 = self.U0
        self.hill_height = self.hill_height
        self.hill_width = self.hill_width
        self.ridge_wavelength = self.ridge_wavelength

        self.R_gas = 287.4
        self.gamm = 1.4
        self.T_ref = 300.0
        self.p_ref = 1e5
        self.h_ref = 10000.0  # [m]
        self.t_ref = 1000.0  # [s]
        self.u_ref = self.h_ref / self.t_ref  # 10 m/s

        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.gravity_direction = 1
        self.gravity_strength = np.zeros(3)
        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.i_gravity = np.zeros(3)
        self.i_gravity[1] = 1

        # x: +-50 km periodic (h(+-50 km) ~ e^-100 h0, seam-safe);
        # y: 0..19.6 km, sponge above ~11.6 km (two vertical wavelengths
        # 2 pi U / N = 6.3 km below it)
        self.xmin = -50000.0 / self.h_ref
        self.xmax = 50000.0 / self.h_ref
        self.ymin = 0.0
        self.ymax = 19600.0 / self.h_ref
        self.zmin = 0.0
        self.zmax = 1.0

        self.u_wind_speed = self.U0 / self.u_ref
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC,
            opts.BdryType.WALL,  # switched to RAYLEIGH in sol_init
            opts.BdryType.PERIODIC,
        )

        self.rayleigh_bdry_switch = True
        self.rayleigh_forcing = False

        self.is_compressible = 1
        self.is_nonhydrostatic = 1

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9
        self.dtfixed = 12.5 / self.t_ref  # U dt / dx ~ 0.32
        self.dtfixed0 = 12.5 / self.t_ref

        self.inx = 256 + 1  # dx ~ 390 m (lambda / dx ~ 10)
        self.iny = 64 + 1  # dy ~ 306 m (lambda_z / dy ~ 21)
        self.inz = 1  # NATIVE 2D (lap2D terrain path)

        # sponge: cells above ~11.6 km
        self.inbcy = 26

        # SLEVE decay scales (Schär et al. 2002): s1 = 15 km, s2 = 2.5 km
        self.vertical_transform = terrain.SLEVETransform(
            s1=15000.0 / self.h_ref, s2=2500.0 / self.h_ref
        )

        # regression run: 300 s spin-up (deterministic gate — exercises the
        # native-2D terrain BC, metric advection, SLEVE elliptic tensor,
        # sponge); the analytic oracle integrates further on a reduced grid
        self.tout = [1e6]
        self.stepmax = 24

        self.stratification = self.stratification_function
        self.orography = self.orography_function
        self.orography_grad = (self.orography_dx, self.orography_zero_grad)
        self.orography_smooth = self.orography_smooth_function
        self.orography_smooth_grad = (
            self.orography_smooth_dx,
            self.orography_zero_grad,
        )
        self.output_timesteps = True

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_schaer_ridge"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)
        self.autogen_fn = False

        self.diag_state = make_diag_state(
            "test_schaer_ridge",
            "target_schaer_ridge",
            self.inx,
            self.iny,
            self.stepmax,
            plot_compare=True,
            # terrain elliptic solves amplify cross-platform rounding: a CI
            # runner deviated 2.4e-5 on rhou from a locally generated target
            # (local same-machine scatter is ~1e-7). Gate well above platform
            # noise; physics is guarded by the linear FFT oracle + SLEVE
            # discriminator.
            tolerances={k: 5e-4 for k in _COMPARED_FIELDS},
        )

    def stratification_function(self, y):
        # constant buoyancy frequency: theta = exp(N^2 z / g)
        Nsq = (self.NN * self.t_ref) ** 2
        g = self.gravity_strength[1] / self.Msq
        return np.exp(Nsq * y / g)

    # --- orography: total, smooth part, analytic gradients (nondim) ---------

    def _scales(self):
        h0 = self.hill_height / self.h_ref
        a = self.hill_width / self.h_ref
        lam = self.ridge_wavelength / self.h_ref
        return h0, a, lam

    def orography_function(self, xi1, xi2):
        h0, a, lam = self._scales()
        return h0 * np.exp(-((xi1 / a) ** 2)) * np.cos(np.pi * xi1 / lam) ** 2 + (
            0.0 * xi2
        )

    def orography_dx(self, xi1, xi2):
        h0, a, lam = self._scales()
        env = np.exp(-((xi1 / a) ** 2))
        return h0 * env * (
            -(2.0 * xi1 / a**2) * np.cos(np.pi * xi1 / lam) ** 2
            - (np.pi / lam) * np.sin(2.0 * np.pi * xi1 / lam)
        ) + (0.0 * xi2)

    def orography_smooth_function(self, xi1, xi2):
        h0, a, _ = self._scales()
        return 0.5 * h0 * np.exp(-((xi1 / a) ** 2)) + 0.0 * xi2

    def orography_smooth_dx(self, xi1, xi2):
        h0, a, _ = self._scales()
        return -h0 * xi1 / a**2 * np.exp(-((xi1 / a) ** 2)) + 0.0 * xi2

    @staticmethod
    def orography_zero_grad(xi1, xi2):
        return 0.0 * xi1 + 0.0 * xi2


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    apply_rayleigh_bdry(ud, elem, node, with_tau=True)

    # constant-N background at physical height via fine-grid quadrature
    hydrostatics.integrated_state(npf, elem, node, th, ud)

    S0c = npf.HydroState.get_S0c(elem)
    rhoY0 = npf.HydroState.rhoY0  # field mode: full per-column fields

    Sol.rhoY[...] = rhoY0
    Sol.rho[...] = rhoY0 * S0c
    Sol.rhou[...] = Sol.rho * ud.u_wind_speed
    Sol.rhov[...] = 0.0
    Sol.rhow[...] = 0.0
    Sol.rhoX[...] = Sol.rho * (Sol.rho / Sol.rhoY - S0c)

    # hydrostatically balanced background: zero perturbation pressure
    npf.p2_nodes[...] = 0.0

    return Sol
