"""3D Agnesi mountain wave — first true-3D terrain golden master (G2 != 0).

The canonical 3D extension of ``test_agnesi_hydrostatic``: an isolated
circular Agnesi bell h(x, z) = h0 / (1 + (x^2 + z^2)/a^2)^(3/2) in uniform
flow U along x, constant buoyancy frequency N = 0.01 1/s, hydrostatic-regime
parameters (N a / U = 10) and linear amplitude (N h0 / U = 0.1). Gravity is
on the default axis (vertical = y), so the orography is a genuine function
of BOTH horizontal coordinates and the second terrain slope G2 is nonzero —
this is the case that locks the G2-carrying metric terms (contravariant
vertical flux, pressure-gradient map, elliptic cross tensor, bottom BC) as
a golden master. Rayleigh sponge above 12 km absorbs the upward-radiating
waves, exactly as in the 2D case.

The regression run (stepmax steps of transient spin-up) is a determinism
gate for the G2 != 0 dynamics; the physics of the shared machinery is
proven by the Smith (1980) oracle on the 2D configuration
(``test_scripts/test_agnesi_analytic.py``). Deliberately NOT wired into CI:
the full-3D elliptic solves make it minutes-scale, too slow for the runner.
"""

import numpy as np

from ..utils import options as opts
from ..flow_solver.physics import hydrostatics

from .case_setup import apply_rayleigh_bdry, build_bdry, make_diag_state

_COMPARED_FIELDS = ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX", "p2_nodes")


class UserData(object):
    grav = 9.81  # [m/s^2]
    omega = 0.0

    NN = 0.01  # [1/s] buoyancy frequency
    U0 = 10.0  # [m/s] background wind
    hill_height = 100.0  # [m]
    hill_width = 10000.0  # [m]

    def __init__(self):
        self.grav = self.grav
        self.omega = self.omega
        # instance-ify the case parameters (vars(UserData()) only carries
        # instance attributes onto UserDataInit)
        self.NN = self.NN
        self.U0 = self.U0
        self.hill_height = self.hill_height
        self.hill_width = self.hill_width

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

        # x, z: +-80 km (= +-8 a) periodic — the 3D bell decays like
        # (r/a)^-3, h(8a) ~ 2e-3 h0 ~ 0.2 m, so the periodic wrap of the
        # orography (and its analytic gradient) is smooth at the seams;
        # y: 0..24 km with the sponge above 12 km (two vertical
        # wavelengths lambda_z = 2 pi U / N ~ 6.3 km below it)
        self.xmin = -80000.0 / self.h_ref
        self.xmax = 80000.0 / self.h_ref
        self.ymin = 0.0
        self.ymax = 24000.0 / self.h_ref
        self.zmin = -80000.0 / self.h_ref
        self.zmax = 80000.0 / self.h_ref

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
        self.dtfixed = 50.0 / self.t_ref  # N dt = 0.5
        self.dtfixed0 = 50.0 / self.t_ref

        self.inx = 64 + 1  # dx = 2.5 km (a / dx = 4)
        self.iny = 32 + 1  # dy = 750 m (lambda_z / dy ~ 8.4)
        self.inz = 64 + 1  # dz = 2.5 km (a / dz = 4)

        # sponge: cells above 12 km
        self.inbcy = 16

        # regression run: 500 s of transient spin-up — deterministic gate
        # on the G2 != 0 path (full-3D terrain elliptic tensor, metric
        # advection in both horizontal directions, terrain bottom BC)
        self.tout = [1e6]
        self.stepmax = 10

        self.stratification = self.stratification_function
        self.orography = self.orography_function
        self.orography_grad = (self.orography_dxi1, self.orography_dxi2)
        self.output_timesteps = True

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_agnesi_3d"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)
        self.autogen_fn = False

        self.diag_state = make_diag_state(
            "test_agnesi_3d",
            "target_agnesi_3d",
            self.inx,
            self.iny,
            self.stepmax,
            # 3D fields are contour-plotted as the transverse mid-slice
            plot_compare=True,
            # same gate as test_agnesi_hydrostatic: terrain elliptic solves
            # (O(100) bicgstab iters/step) amplify cross-platform rounding
            # well above the 1e-5 default; local same-machine scatter is
            # ~1e-7. Physics is guarded by the Smith (1980) oracle on the
            # shared 2D configuration.
            tolerances={k: 5e-4 for k in _COMPARED_FIELDS},
        )

    def stratification_function(self, y):
        # constant buoyancy frequency: theta = exp(N^2 z / g)
        Nsq = (self.NN * self.t_ref) ** 2
        g = self.gravity_strength[1] / self.Msq
        return np.exp(Nsq * y / g)

    def orography_function(self, xi1, xi2):
        # circular 3D Agnesi bell: h0 / (1 + r^2/a^2)^(3/2)
        h0 = self.hill_height / self.h_ref
        a = self.hill_width / self.h_ref
        return h0 / (1.0 + (xi1**2 + xi2**2) / a**2) ** 1.5

    def orography_dxi1(self, xi1, xi2):
        h0 = self.hill_height / self.h_ref
        a = self.hill_width / self.h_ref
        return -3.0 * h0 * xi1 / a**2 / (1.0 + (xi1**2 + xi2**2) / a**2) ** 2.5

    def orography_dxi2(self, xi1, xi2):
        h0 = self.hill_height / self.h_ref
        a = self.hill_width / self.h_ref
        return -3.0 * h0 * xi2 / a**2 / (1.0 + (xi1**2 + xi2**2) / a**2) ** 2.5


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
