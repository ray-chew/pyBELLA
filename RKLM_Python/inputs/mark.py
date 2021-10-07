import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
from physics.hydrostatics import hydrostatic_state, hydrostatic_column, hydrostatic_initial_pressure
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2
from management.variable import States


class UserData(object):
    NSPEC = 1

    grav = 9.81                 # [m s^{-2}]
    omega = 7.292 * 1e-5        # [s^{-1}]
    omega = 0.0

    R_gas = 287.4
    R_vap = 461.0
    Q_vap = 2.53e+06
    gamma = 1.4
    cp_gas = gamma * R_gas / (gamma-1.0)

    p_ref = 1e+5
    T_ref = 300.00              # [K]
    rho_ref = p_ref / (R_gas * T_ref)
    N_ref = grav / np.sqrt(cp_gas * T_ref)
    Cs = np.sqrt(gamma * R_gas * T_ref)

    # h_ref = 20.0e3
    h_ref = 10.0e3              # [m]
    # t_ref = h_ref / Cs
    t_ref = 100.0               # [s]
    u_ref = h_ref / t_ref

    i_gravity = np.zeros((3))
    i_coriolis = np.zeros((3))

    def __init__(self):
        self.h_ref = self.h_ref
        self.t_ref = self.t_ref
        self.T_ref = self.T_ref
        self.p_ref = self.p_ref
        self.rho_ref = self.rho_ref
        self.u_ref = self.u_ref
        self.Nsq_ref = self.N_ref * self.N_ref
        self.g_ref = self.grav
        self.gamm = self.gamma
        self.Rg_over_Rv = self.R_gas / self.R_vap
        self.Q = self.Q_vap / (self.R_gas * self.T_ref)
        self.Rg = self.R_gas / (self.h_ref**2 / self.t_ref**2 / self.T_ref)
        self.cp_gas = self.cp_gas

        self.nspec = self.NSPEC

        self.is_nonhydrostatic = 1
        self.is_compressible = 1
        self.is_ArakawaKonor = 0

        self.compressibility = 1.0
        self.acoustic_timestep = 0
        self.acoustic_order = 0
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.gravity_strength = np.zeros((3))
        self.coriolis_strength = np.zeros((3))

        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.coriolis_strength[1] = self.omega * self.t_ref

        for i in range(3):
            if (self.gravity_strength[i] > np.finfo(np.float).eps) or (i == 1):
                self.i_gravity[i] = 1
                self.gravity_direction = 1

            if (self.coriolis_strength[i] > np.finfo(np.float).eps):
                self.i_coriolis[i] = 1

        # self.xmin = - 24.0e6 / self.h_ref
        # self.xmax =   24.0e6 / self.h_ref
        j = 2.0
        Lx = 1.0 * np.pi * self.Cs / self.N_ref * j
        self.xmin = - Lx / self.h_ref
        self.xmax =   Lx / self.h_ref
        self.ymin = - 0.0
        self.ymax =   4.0 * 8776.0 / self.h_ref
        # self.ymax =   1.0
        self.zmin = - 1.0
        self.zmax =   1.0

        self.u_wind_speed = 0.0 * self.u_ref
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.WALL
        self.bdry_type[2] = BdryType.WALL

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9/2.0
        # self.CFL = 0.9

        self.inx = 301+1
        self.iny = 10+1
        self.inz = 1

        self.dtfixed0 = 0.5 * ((self.xmax - self.xmin) / (self.inx-1)) / 1.0
        self.dtfixed = self.dtfixed0

        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.tol = 1.e-8
        self.max_iterations = 6000

        self.perturb_type = 'pos_perturb'
        self.blending_mean = 'rhoY' # 1.0, rhoY
        self.blending_conv = 'rho' #theta, rho
        self.blending_type = 'half' # half, full

        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.blending_weight = 0./16

        self.initial_blending = False

        self.initial_projection = True
        self.initial_impl_Euler = False

        # self.tout = np.arange(0.0,9000.0,1000.0)[1:]
        # self.tout = np.arange(0.0,205.0,5.0)[1:]
        self.tout = [8000.0]
        self.stepmax = 20000

        self.output_base_name = "_mark_wave"

        aux = 'debug_ic'
        self.aux = aux

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_function
        self.output_timesteps = True

    def stratification_function(self, y):
        Nsq = self.Nsq_ref * self.t_ref**2
        g = self.gravity_strength[1] / self.Msq

        return np.exp(Nsq * y / g)

    def rhoe_function(self,rho,u,v,w,p,ud,th):
        pass


def sol_init(Sol, mpv, elem, node, th, ud, seeds=None):
    def bump(xi):
        # eqn (11)
        tmp = np.zeros_like(xi)
        tmp[np.where(np.abs(xi) < 1.0)] = np.exp(-1.0 / (1.0 - xi[np.where(np.abs(xi) < 1.0)]**2))
        # return tmp
        return np.ones_like(xi)

    hydrostatic_state(mpv, elem, node, th, ud)
    n = 4.0           # eqn (12)
    A0 = 1.0e-3     # eqn (12)

    # ud.Msq = np.sqrt(ud.Msq)

    Msq = ud.Msq
    g = ud.gravity_strength[1]
    kappa = th.Gamma
    Omega = ud.coriolis_strength[1]

    waveno = n * 2.0 * np.pi / (0.5 * (ud.xmax - ud.xmin))

    x = elem.x.reshape(-1,1)
    y = elem.y.reshape(1,-1)


    Hrho = 1.0 / g
    kGam = (1.0 - th.gamm / 2.0) / Hrho                 # eqn (8)
    A = A0 * bump(2.0 * y / (n * Hrho) - 1)             # eqn (12)

    use_hydrostate = False
    
    if use_hydrostate:
        # Use hydrostatically balanaced background
        rhobar = mpv.HydroState.rho0.reshape(1,-1)
        Ybar = mpv.HydroState.Y0.reshape(1,-1)
        pibar = mpv.HydroState.p20.reshape(1,-1)
        div = 1.0
    else:
        # Use hydrostatic balance in Mark's notes
        Htheta = Hrho / kappa
        rhobar = np.exp(-y / Hrho) # eqn (7)
        Ybar = np.exp(y / Htheta)  # eqn (9)
        pibar = 1.0 / Ybar
        div = Msq

    N = ud.t_ref * np.sqrt(ud.Nsq_ref)                  # eqn (5) dimensionless Brunt-Väisälä frequency
    Cs = np.sqrt(th.gamm / Msq)                         # eqn (6) dimensionless speed of sound

    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed

    # 1.0 / c_p in the notes
    # fac = 1.0 / (ud.cp_gas / (ud.h_ref**2 / ud.t_ref**2 / ud.T_ref))
    fac = th.Gamma


    # set up perturbation quantities
    # exp(-kGam * y)  / sqrt(rhobar) / Ybar = 1.0
    up = A * np.sqrt(1.0 / rhobar) * np.exp(-kGam * y) * np.cos(N / Cs * waveno * x)
    vp = -A * np.sqrt(Omega / (Cs * kGam)) * np.sqrt(1.0 / rhobar) * np.exp(-kGam * y) * np.sin(N / Cs * waveno * x)
    wp = 0.0
    Yp = A * np.sqrt(Omega / (Cs * kGam)) * N / g * Ybar * np.sqrt(1.0 / rhobar) * np.exp(-kGam * y) * np.sin(N / Cs * waveno * x)     # eqn (3)
    pi = A * Cs * fac * np.sqrt(1.0 / rhobar) / Ybar * np.exp(-kGam * y) * np.cos(N / Cs * waveno * x)     # eqn (4)

    u = u0 + up
    v = v0 + vp
    w = w0 + wp
    Y = Ybar + Yp

    rho = ((pibar * Msq / div + pi)**th.gm1inv) / Y

    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhow[...] = rho * w
    Sol.rhoe[...] = 0.0
    Sol.rhoY[...] = rho * Y
    Sol.rhoX[...] = 0.0
    mpv.p2_cells[...] = pi / Msq

    ###################################################
    # initialise nodal pi
    xn = node.x.reshape(-1,1)
    yn = node.y.reshape(1,-1)

    # initialise nodal pressure
    Hrho_n = Hrho
    kGam_n = (1.0 - th.gamm / 2.0) / Hrho_n

    if use_hydrostate:
        # Use hydrostatically balanced background
        rhobar_n = mpv.HydroState_n.rho0.reshape(1,-1)
        Ybar_n = mpv.HydroState_n.Y0.reshape(1,-1)
        pibar_n = mpv.HydroState_n.p20.reshape(1,-1)
    else:
        # Use hydrostatic balance from notes
        rhobar_n = np.exp(-yn / Hrho_n)
        Ybar_n = np.exp(-yn / Htheta)
        pibar_n = 1.0 / Ybar_n

    An = A0 * bump(2.0 * yn / (n * Hrho_n) - 1.0)

    pi_n = An * Cs * fac * np.sqrt(1.0 / rhobar_n) / Ybar_n * np.exp(-kGam_n * yn) * np.cos(N / Cs * waveno * xn)

    mpv.p2_nodes[...] = pi_n / Msq

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    return Sol
