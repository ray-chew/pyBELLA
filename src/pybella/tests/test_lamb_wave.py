import numpy as np

from ..utils import options as opts
from ..utils.data_structures import DiagnosticState

from ..flow_solver.utils import boundary as bdry
from ..flow_solver.physics import hydrostatics


class UserData(object):

    def __init__(self):
        self.grav = 9.81
        self.h_ref = 10.0e3  # [m]
        self.t_ref = 100.0  # [s]
        self.T_ref = 300.00  # [K]
        self.p_ref = 1e5
        self.u_ref = self.h_ref / self.t_ref
        self.R_gas = 287.4
        self.gamm = 1.4
        self.Cs = np.sqrt(self.gamm * self.R_gas * self.T_ref)
        self.cp_gas = self.gamm * self.R_gas / (self.gamm - 1.0)
        self.N_ref = self.grav / np.sqrt(self.cp_gas * self.T_ref)
        self.Rg = self.R_gas / (self.h_ref**2 / self.t_ref**2 / self.T_ref)

        self.is_nonhydrostatic = 1
        self.is_compressible = 1
        self.is_ArakawaKonor = 0

        self.compressibility = 1.0
        self.acoustic_timestep = 0
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.gravity_strength = np.zeros((3))
        self.coriolis_strength = np.zeros((3))

        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)

        gravity_mask = (self.gravity_strength > np.finfo(np.float64).eps) | (
            np.arange(3) == 1
        )
        self.i_gravity = gravity_mask.astype(int)
        if np.any(gravity_mask):
            self.gravity_direction = np.where(gravity_mask)[0][
                -1
            ]  # Use last matching index

        j = 4.0
        Lx = 1.0 * np.pi * self.Cs / self.N_ref * j
        self.xmin = -Lx / self.h_ref
        self.xmax = Lx / self.h_ref
        self.ymin = -0.0
        self.ymax = 8.0
        self.zmin = -1.0
        self.zmax = 1.0

        self.ATMOSPHERIC_EXTENSION = True
        self.rayleigh_bdry_switch = False

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9

        self.inx = 151 + 1
        self.iny = 60 + 1
        self.inz = 1

        self.dtfixed0 = 10.0 / self.t_ref
        self.dtfixed = self.dtfixed0

        self.tol = 1.0e-30
        self.max_iterations = 10000

        self.tout = [360.0]
        # self.tout = np.arange(0,361,1.0)
        # self.tout = np.append(self.tout, [720.0])
        self.stepmax = 301
        self.output_timesteps = True

        self.autogen_fn = False

        self.diag = True
        self.diag_updt_targets = False

        self.output_base_name = "_lamb_wave"
        self.output_type = "test" if not self.diag_updt_targets else "target"
        self.aux = ""
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.iny - 1)

        self.diag_state = DiagnosticState(
            test_name="test_lamb_wave",
            file_name="target_lamb_wave",
            Nx=self.inx - 1,
            Ny=self.iny - 1,
            steps=[self.stepmax - 1],
        )

        self.stratification = self.stratification_wrapper
        self.init_forcing = self.forcing

    def stratification_wrapper(self, dy):
        return lambda y: self.stratification_function(y, dy)

    def stratification_function(self, y, dy):
        g = self.gravity_strength[1]

        Gamma = (self.gamm - 1.0) / self.gamm
        Hex = 1.0 / (Gamma * g)
        pi_m = np.exp(-(y - 0.5 * dy) / Hex)
        pi_p = np.exp(-(y + 0.5 * dy) / Hex)

        Theta = -(Gamma * g * dy) / (pi_p - pi_m)

        return Theta

    @staticmethod
    def rayleigh_bc_function(ud):
        if ud.bdry_type[1] == opts.BdryType.RAYLEIGH or ud.rayleigh_forcing == True:
            ud.inbcy = ud.iny - 1
            ud.iny0 = np.copy(ud.iny)
            ud.iny = ud.iny0 + int(3 * ud.inbcy)

            # tentative workaround
            ud.bcy = ud.ymax
            ud.ymax += 3.0 * ud.bcy

    class forcing(object):
        def __init__(
            self,
            k,
            mu,
            Cs,
            F,
            N,
            Gamma,
            ampl,
            g,
            rhobar,
            Ybar,
            rhobar_n,
            Ybar_n,
            X,
            Y,
            Xn,
            Yn,
        ):
            self.k = k
            self.mu = mu
            self.Cs = Cs
            self.F = F
            self.N = N
            self.Gamma = Gamma

            self.oorhobarsqrt = 1.0 / np.sqrt(rhobar.T)
            self.Ybar = Ybar.T

            self.oorhobarsqrt_n = 1.0 / np.sqrt(rhobar_n.T)
            self.Ybar_n = Ybar_n.T

            self.g = g
            self.ampl = ampl

            self.X = X
            self.Y = Y
            self.Xn = Xn
            self.Yn = Yn

        def get_T_matrix(self):
            # system matrix of linearized equations
            matrix = -np.array(
                [
                    [0, self.F, 0, 1j * self.Cs * self.k],
                    [-self.F, 0, -self.N, self.Cs * (self.mu + self.Gamma)],
                    [0, self.N, 0, 0],
                    [1j * self.Cs * self.k, self.Cs * (self.mu - self.Gamma), 0, 0],
                ]
            )

            self.T_matrix = matrix

        def eigenfunction(self, t, s, grid="c"):
            if grid == "c":
                x, z = self.X, self.Y
            elif grid == "n":
                (
                    x,
                    z,
                ) = (
                    self.Xn,
                    self.Yn,
                )

            # Compute eigenvalues and eigenvectors
            eigval, eigvec = np.linalg.eig(self.T_matrix)

            # Find index of eigenvalue
            # with greatest real part aka the instability growth rate
            ind = np.argmax(np.real(eigval))

            # construct solution according to eq. 2.27 and 2.19
            exponentials = np.exp(
                1j * self.k * x + self.mu * z + (eigval[ind]) * (t) + 1j * s * t
            )
            chi_u = self.ampl * np.real(eigvec[0, ind] * exponentials)
            chi_w = self.ampl * np.real(eigvec[1, ind] * exponentials)
            chi_th = self.ampl * np.real(eigvec[2, ind] * exponentials)
            chi_pi = self.ampl * np.real(eigvec[3, ind] * exponentials)

            self.arrs = (chi_u, chi_w, chi_th, chi_pi)

        def dehatter(self, th, grid="c"):
            if grid == "n":
                Ybar = self.Ybar_n
                oorhobarsqrt = self.oorhobarsqrt_n
            elif grid == "c":
                Ybar = self.Ybar
                oorhobarsqrt = self.oorhobarsqrt

            chi_u, chi_v, chi_Y, chi_pi = self.arrs

            up = oorhobarsqrt * chi_u
            vp = oorhobarsqrt * chi_v
            Yp = oorhobarsqrt * self.N / self.g * Ybar * chi_Y
            pi_p = oorhobarsqrt * self.Cs / Ybar / th.Gammainv * chi_pi

            return up.T, vp.T, Yp.T, pi_p.T


def sol_init(Sol, npf, elem, node, th, ud, seeds=None):
    if hasattr(ud, "rayleigh_bdry_switch"):
        if ud.rayleigh_bdry_switch:
            ud.bdry_type[1] = opts.BdryType.RAYLEIGH

    if ud.bdry_type[1] == opts.BdryType.RAYLEIGH:
        ud.tcy, ud.tny = bdry.get_tau_y(ud, elem, node, 0.5)

    A0 = 1.0e-1 / ud.u_ref
    Msq = ud.Msq
    g = ud.gravity_strength[1] * ud.Rg

    x = elem.x.reshape(-1, 1)
    y = elem.y.reshape(1, -1)
    X, Y = np.meshgrid(x, y)

    xn = node.x.reshape(-1, 1)
    yn = node.y.reshape(1, -1)

    dy = np.diff(node.y)[0]

    Xn, Yn = np.meshgrid(xn, yn)

    ##################################################
    # Following Rupert's fix, reinitialise all background quantities
    # as derived from one quantity.
    ud.stratification = ud.stratification(dy)

    # Use hydrostatically balanced background
    hydrostatics.analytical_state(npf, elem, node, th, ud)
    rhobar = npf.HydroState.rho0.reshape(1, -1)
    Ybar = npf.HydroState.Y0.reshape(1, -1)
    pibar = npf.HydroState.p20.reshape(1, -1) * ud.Msq

    rhobar_n = npf.HydroState_n.rho0.reshape(1, -1)
    Ybar_n = npf.HydroState_n.Y0.reshape(1, -1)

    ##################################################
    # dimensionless Brunt-Väisälä frequency
    N = ud.t_ref * np.sqrt(ud.Nsq_ref)
    # dimensionless speed of sound
    Cs = np.sqrt(th.gamm / Msq)
    ud.Cs = Cs
    ud.Ns = N
    # dimensionless Coriolis strength
    if ud.coriolis_strength[2] == 0.0:
        ud.coriolis_strength[2] += 1e-15
    F = ud.coriolis_strength[2]

    G = np.sqrt(9.0 / 40.0)
    Gamma = G * N / Cs
    k = N / Cs

    ud.rf_bot = ud.init_forcing(
        k, -Gamma, Cs, F, N, Gamma, A0, g, rhobar, Ybar, rhobar_n, Ybar_n, X, Y, Xn, Yn
    )
    ud.rf_bot.get_T_matrix()

    ud.u_wind_speed = 0.0

    ud.rf_bot.eigenfunction(0, 1)
    up, vp, Yp, pi_p = ud.rf_bot.dehatter(th)

    u = ud.u_wind_speed + up
    v = ud.v_wind_speed + vp
    w = ud.w_wind_speed
    Y = Ybar + Yp

    rho = rhobar

    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhow[...] = rho * w
    Sol.rhoY[...] = rho * Y
    Sol.rhoX[...] = 0.0
    npf.p2_cells[...] = pi_p

    ###################################################
    # initialise nodal pi
    ud.rf_bot.eigenfunction(0, 1, grid="n")
    _, _, _, pi_n = ud.rf_bot.dehatter(th, grid="n")

    npf.p2_nodes[...] = pi_n

    bdry.set_explicit_boundary_data(Sol, elem, ud, th, npf)

    if hasattr(ud, "mixed_run"):
        if ud.mixed_run:
            ud.coriolis_strength[2] = 2.0 * 7.292 * 1e-5 * ud.t_ref

    if hasattr(ud, "trad_forcing"):
        if ud.trad_forcing:
            ud.rf_bot.F = 0.0
            ud.rf_bot.get_T_matrix()

    return Sol
