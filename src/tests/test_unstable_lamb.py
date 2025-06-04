"""
Unstable Lamb Wave integral test involving vertical Coriolis and Rayleigh BC.
"""

import numpy as np
from ..flow_solver.physics import hydrostatics
from ..flow_solver.utils import boundary as bdry

from ..utils.data_structures import DiagnosticState


class UserData(object):
    def __init__(self):
        self.grav = 9.81  # [m/s^2]
        self.omega = 7.292 * 1e-5  # [s^{-1}]
        self.t_ref = 100.0  # [s]

        self.h_ref = 10.0e3              # [m]
        self.u_ref = self.h_ref / self.t_ref

        self.T_ref = 300.00              # [K]
        self.R_gas = 287.4               # [J kg^{-1} K^{-1}]
        self.gamma = 1.4
        self.cp_gas = self.gamma * self.R_gas / (self.gamma-1.0)

        self.Nsq = (self.grav / np.sqrt(self.cp_gas * self.T_ref))**2
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        ##########################################
        # SPATIAL GRID
        ##########################################
        self.inx = 301 + 1
        self.iny = 30 + 1
        self.inz = 1

        self._setup_domain()

        ##########################################
        # BOUNDARY CONDITIONS
        ##########################################
        from ..flow_solver.utils import options as opts

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = opts.BdryType.PERIODIC
        self.bdry_type[1] = opts.BdryType.WALL
        self.bdry_type[2] = opts.BdryType.WALL
        self.ATMOSPHERIC_EXTENSION = True

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9
        self.dtfixed0 = 100.0 / self.t_ref
        self.dtfixed = self.dtfixed0

        self.tout = [36.0]
        self.stepmax = 101

        self.is_compressible = 1
        self.is_nonhydrostatic = 1
        self.is_ArakawaKonor = 0

        self.compressibility = 1.0
        self.acoustic_timestep = 0

        ##########################################
        # PHYSICS AND BACKGROUND WIND
        ##########################################
        self.u_wind_speed = 0.0
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        ##########################################
        # BLENDING
        ##########################################
        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.blending_weight = 0.0 / 16
        self.blending_mean = "rhoY"
        self.blending_conv = "rho"
        self.blending_type = "half"
        self.initial_blending = False

        ##########################################
        # STRATIFICATION
        ##########################################
        self.stratification = self.stratification_wrapper

        ##########################################
        # DIAGNOSTICS
        ##########################################
        self.diag = True
        self.diag_updt_targets = True

        ##########################################
        # OUTPUTS
        ##########################################
        self.output_base_name = "_unstable_lamb"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.output_timesteps = True

        # Rayleigh forcing parameters
        self.rayleigh_forcing = True
        self.rayleigh_forcing_type = "func"
        self.rayleigh_forcing_fn = None
        self.rayleigh_forcing_path = None

        self.diag_state = DiagnosticState(
            test_name="test_unstable_lamb",
            file_name="target_unstable_lamb",
            Nx=self.inx - 1,
            Ny=self.iny - 1,
            steps=[self.stepmax - 1],
            plot_compare=True,
        )

        self.autogen_fn = False

    def _setup_domain(self):
        """Setup domain dimensions based on wave parameters"""
        # Reference values (from global constants)
        R_gas = 287.4  # [J kg^{-1} K^{-1}]
        gamma = 1.4
        T_ref = 300.0  # [K]
        h_ref = 10.0e3  # [m]

        # Compute derived quantities
        Cs = np.sqrt(gamma * R_gas * T_ref)
        cp_gas = gamma * R_gas / (gamma - 1.0)
        N_ref = self.grav / np.sqrt(cp_gas * T_ref)

        # Domain size calculation
        j = 4.0
        Lx = 1.0 * np.pi * Cs / N_ref * j

        self.xmin = -Lx / h_ref
        self.xmax = Lx / h_ref
        self.ymin = 0.0
        self.ymax = 2.0
        self.zmin = -1.0
        self.zmax = 1.0

    def stratification_wrapper(self, dy):
        return lambda y : self.stratification_function(y, dy)

    def stratification_function(self, y, dy):

        Nsq = self.Nsq * self.t_ref**2
        g = self.grav / self.Msq

        return np.exp(Nsq * y / g)



def sol_init(Sol, mpv, elem, node, th, ud, seeds=None):
    def bump(xi):
        # eqn (11)
        tmp = np.zeros_like(xi)
        tmp[np.where(np.abs(xi) < 1.0)] = np.exp(
            -1.0 / (1.0 - xi[np.where(np.abs(xi) < 1.0)] ** 2)
        )
        # return tmp
        return np.ones_like(xi)

    if ud.bdry_type[1].value == "radiation":
        ud.tcy, ud.tny = bdry.get_tau_y(ud, elem, node, 0.5)

        if hasattr(ud, "rayleigh_forcing"):
            if ud.rayleigh_forcing:
                ud.forcing_tcy, ud.forcing_tny = bdry.get_bottom_tau_y(
                    ud, elem, node, 0.2, cutoff=0.1
                )

    A0 = 1.0e-3  # eqn (12)

    Msq = ud.Msq
    g = ud.gravity_strength[1]
    kappa = th.Gamma

    x = elem.x.reshape(-1, 1)
    y = elem.y.reshape(1, -1)
    dy = np.diff(node.y)[0]

    Hrho = 1.0 / g
    use_hydrostate = False

    ud.stratification = ud.stratification(dy)
    hydrostatics.state(mpv, elem, node, th, ud)

    if use_hydrostate:
        # Use hydrostatically balanaced background
        rhobar = mpv.HydroState.rho0.reshape(1, -1)
        Ybar = mpv.HydroState.Y0.reshape(1, -1)
        pibar = mpv.HydroState.p20.reshape(1, -1) * ud.Msq

    else:
        # Use hydrostatic balance in Mark's notes
        Htheta = Hrho / kappa
        rhobar = np.exp(-y / Hrho)
        Ybar = np.exp(y / Htheta)
        pibar = 1.0 / Ybar
        mpv.HydroState.rho0[...] = rhobar
        mpv.HydroState.Y0[...] = Ybar
        mpv.HydroState.S0[...] = 1.0 / Ybar  # 1.0 / ud.stratification(y)
        mpv.HydroState.rhoY0[...] = rhobar * Ybar

    # dimensionless Brunt-Väisälä frequency
    N = ud.t_ref * np.sqrt(ud.Nsq_ref)
    # dimensionless speed of sound
    Cs = np.sqrt(th.gamm / Msq)
    waveno = N / Cs

    ud.u_wind_speed = 0.0  # -Cs

    Gamma = 1.0 / Hrho * (1.0 / th.gamm - 0.5)

    # time shift of the initial solution
    ts = -0.5 / N * np.pi
    ts = 0.0

    # set up perturbation quantities
    # exp(-kGam * y)  / sqrt(rhobar) / Ybar = 1.0
    up = A0 * Ybar * np.cos(waveno * x + Cs * ts)
    vp = 0.0
    wp = 0.0
    Yp = 0.0
    # th.Gamma * ud.Msq == 1 / dimensionless(c_p)
    fac = th.Gamma
    pi = A0 * Cs * fac * np.cos(waveno * x + Cs * ts) * Msq

    # up = A / rhobar**0.5 * np.exp(-Gamma * y) * np.cos(N / Cs * (waveno * x + Cs * ts))
    # pi = A * Cs * fac / rhobar**0.5 / Ybar * np.exp(-Gamma * y) * np.cos(N / Cs * (waveno * x + Cs * ts))

    u = ud.u_wind_speed + up
    v = ud.v_wind_speed + vp
    w = ud.w_wind_speed + wp
    Y = Ybar + Yp
    dPdpi = th.gm1inv * pibar ** (th.gm1inv - 1.0)
    Pbar = pibar**th.gm1inv
    rhoY = Pbar  # + dPdpi * pi

    # eqn (2.3)
    # rho = (((pibar + pi))**th.gm1inv) / Y
    # rho = Pbar / Y
    rhobar = Pbar / Ybar
    rho = rhobar  # rhoY / Ybar

    Sol.rho[...] = Pbar / Y
    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhow[...] = rho * w
    Sol.rhoY[...] = rhoY
    Sol.rhoX[...] = 0.0
    mpv.p2_cells[...] = pi / Msq

    ###################################################
    # initialise nodal pi
    xn = node.x.reshape(-1, 1)
    yn = node.y.reshape(1, -1)

    # initialise nodal pressure
    Hrho_n = Hrho

    if use_hydrostate:
        # Use hydrostatically balanced background
        Ybar_n = mpv.HydroState_n.Y0.reshape(1, -1)
        rhobar_n = mpv.HydroState_n.rho0.reshape(1, -1)
    else:
        # Use hydrostatic balance from notes
        Ybar_n = np.exp(yn / Htheta)
        rhobar_n = np.exp(-yn / Hrho_n)
        mpv.HydroState_n.Y0[...] = Ybar_n
        mpv.HydroState_n.S0[...] = 1.0 / Ybar_n

    An = A0

    pi_n = An * Cs * fac * np.cos(waveno * xn + Cs * ts)
    # pi_n = An * Cs * fac / rhobar_n**0.5 / Ybar_n * np.exp(-Gamma * yn) * np.cos(N / Cs * (waveno * xn + Cs * ts))

    mpv.p2_nodes[...] = pi_n

    # if ud.bdry_type[1] == 'RAYLEIGH':
    #     rayleigh_damping(Sol, mpv, ud, ud.tcy, elem, th)

    rhoY_tmp = np.copy(Sol.rhoY)
    rho_tmp = np.copy(Sol.rho)

    bdry.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    return Sol
