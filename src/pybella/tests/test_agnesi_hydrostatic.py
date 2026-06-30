"""Agnesi hydrostatic mountain wave — terrain-following golden-master case.

Linear hydrostatic regime (Smith 1980): constant buoyancy frequency
N = 0.01 1/s, uniform wind U = 10 m/s over a witch-of-Agnesi hill with
half-width a = 10 km (N a / U = 10, hydrostatic) and height h0 = 100 m
(N h0 / U = 0.1, linear). Gal-Chen terrain-following coordinates, quasi-2D
3D vertical slice (z degenerate periodic) through the full-tensor elliptic
operator. Rayleigh sponge above 12 km (two vertical wavelengths
lambda_z = 2 pi U / N ~ 6.3 km below it) absorbs upward-radiating waves.

The regression run (stepmax steps) is a determinism gate; the physics is
proven by ``test_scripts/test_agnesi_analytic.py``, which integrates to a
quasi-steady state and compares against the Smith (1980) analytic wave
field and wave drag (``tests/agnesi_smith_analytic.py``).
"""

import numpy as np

from ..utils import options as opts
from ..utils import axes
from ..flow_solver.physics import hydrostatics
from ..flow_solver.utils.boundary import rayleigh_boundary as bdry_r

from .case_setup import build_bdry, make_diag_state

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
        # instance-ify the case parameters: vars(UserData()) only carries
        # instance attributes onto UserDataInit, and the Smith comparator
        # reads them from ud
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

        # x: +-100 km (= +-10 a) periodic — hydrostatic waves radiate
        # vertically, so the wave field stays near the hill and the top
        # sponge absorbs it; y: 0..24 km with the sponge above 12 km
        self.xmin = -100000.0 / self.h_ref
        self.xmax = 100000.0 / self.h_ref
        self.ymin = 0.0
        self.ymax = 24000.0 / self.h_ref
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
        self.dtfixed = 50.0 / self.t_ref  # N dt = 0.5
        self.dtfixed0 = 50.0 / self.t_ref

        self.inx = 128 + 1  # dx = 1.5625 km (a / dx = 6.4)
        self.iny = 64 + 1  # dy = 375 m (lambda_z / dy ~ 17)
        self.inz = 1 + 1

        # sponge: cells above 12 km
        self.inbcy = 32

        # regression run: 750 s (transient spin-up phase, deterministic —
        # exercises terrain BC, metric advection, elliptic tensor, sponge);
        # the analytic oracle extends grid/stepmax to reach steady state
        self.tout = [1e6]
        self.stepmax = 15

        self.stratification = self.stratification_function
        self.orography = self.orography_function
        self.output_timesteps = True

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_agnesi_hydrostatic"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)
        self.autogen_fn = False

        self.diag_state = make_diag_state(
            "test_agnesi_hydrostatic",
            "target_agnesi_hydrostatic",
            self.inx,
            self.iny,
            self.stepmax,
            plot_compare=True,
            # terrain elliptic solves (~130 bicgstab iters/step) amplify
            # cross-platform rounding: first CI run (2026-06-11) deviated
            # 7.0e-5 on rhou / 5.1e-6 on rho from the locally generated
            # target (local same-machine scatter is ~1e-7). Gate well above
            # platform noise; physics is guarded by the Smith (1980) oracle.
            tolerances={k: 5e-4 for k in _COMPARED_FIELDS},
        )

    def stratification_function(self, y):
        # constant buoyancy frequency: theta = exp(N^2 z / g)
        Nsq = (self.NN * self.t_ref) ** 2
        g = self.gravity_strength[1] / self.Msq
        return np.exp(Nsq * y / g)

    def orography_function(self, xi1, xi2):
        h0 = self.hill_height / self.h_ref
        a = self.hill_width / self.h_ref
        return h0 * a**2 / (xi1**2 + a**2) + 0.0 * xi2


def sol_init(Sol, npf, elem, node, th, ud, seed=None):
    if getattr(ud, "rayleigh_bdry_switch", False):
        ud.bdry_type[axes.vertical_axis(ud)] = opts.BdryType.RAYLEIGH
        ud.tcy, ud.tny = bdry_r.get_tau_y(ud, elem, node, 0.5)

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
