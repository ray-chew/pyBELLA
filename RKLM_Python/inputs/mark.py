import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import LimiterType
from physics.hydrostatics import hydrostatic_state
from inputs.boundary import set_explicit_boundary_data, get_tau_y, get_bottom_tau_y, rayleigh_damping

class UserData(object):
    NSPEC = 1
    grav = 9.81                 # [m s^{-2}]
    omega = 7.292 * 1e-5        # [s^{-1}]

    R_gas = 287.4               # [J kg^{-1} K^{-1}]
    R_vap = 461.0
    Q_vap = 2.53e+06
    gamma = 1.4
    cp_gas = gamma * R_gas / (gamma-1.0)

    p_ref = 1e+5
    T_ref = 300.00              # [K]
    rho_ref = p_ref / (R_gas * T_ref)
    N_ref = grav / np.sqrt(cp_gas * T_ref)
    Cs = np.sqrt(gamma * R_gas * T_ref)

    h_ref = 10.0e3              # [m]
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
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.gravity_strength = np.zeros((3))
        self.coriolis_strength = np.zeros((3))

        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.coriolis_strength[2] = 2.0 * self.omega * self.t_ref

        for i in range(3):
            if (self.gravity_strength[i] > np.finfo(np.float).eps) or (i == 1):
                self.i_gravity[i] = 1
                self.gravity_direction = 1

            if (self.coriolis_strength[i] > np.finfo(np.float).eps):
                self.i_coriolis[i] = 1

        j = 4.0
        Lx = 1.0 * np.pi * self.Cs / self.N_ref * j
        self.xmin = - Lx / self.h_ref
        self.xmax =   Lx / self.h_ref
        self.ymin = - 0.0
        self.ymax =   2.5
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
        self.CFL = 0.9

        self.inx = 301+1
        self.iny = 30+1
        self.inz = 1

        self.dtfixed0 = 600.0 / self.t_ref
        self.dtfixed = self.dtfixed0
        
        self.CFL = 0.45
        self.dtfixed0 = 0.5*((self.xmax-self.xmin)/(self.inx-1))/1.0
        self.dtfixed = self.dtfixed0

        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.tol = 1.e-16
        self.max_iterations = 1000

        # blending parameters
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


        self.tout = [720.0]
        self.stepmax = 31
        self.output_timesteps = True

        self.output_base_name = "_mark_wave"
        aux = 'bdl_test_S600_a1'
        self.aux = aux
        self.output_suffix = '_301_120_720.000000_rstrt_init_S600_a1'

        self.stratification = self.stratification_function
        self.rayleigh_bc = self.rayleigh_bc_function

        self.rayleigh_forcing = False
        self.rayleigh_forcing_fn = 'output_mark_wave_ensemble=1_301_120_bottom_forcing_S400.h5'
        self.rayleigh_forcing_path = './output_mark_wave'
        

        # self.rayleigh_bc(self)

    def stratification_function(self, y):
        Nsq = self.Nsq_ref * self.t_ref**2
        g = self.gravity_strength[1] / self.Msq

        return np.exp(Nsq * y / g)

    @staticmethod
    def rayleigh_bc_function(ud):
        if ud.bdry_type[1] == BdryType.RAYLEIGH:
            ud.inbcy = ud.iny - 1
            ud.iny0 = np.copy(ud.iny)
            ud.iny = ud.iny0 + int(3*ud.inbcy)

            # tentative workaround
            ud.bcy = ud.ymax
            ud.ymax += 3.0 * ud.bcy


def sol_init(Sol, mpv, elem, node, th, ud, seeds=None):
    def bump(xi):
        # eqn (11)
        tmp = np.zeros_like(xi)
        tmp[np.where(np.abs(xi) < 1.0)] = np.exp(-1.0 / (1.0 - xi[np.where(np.abs(xi) < 1.0)]**2))
        # return tmp
        return np.ones_like(xi)


    if ud.bdry_type[1].value == 'radiation':
        ud.tcy, ud.tny = get_tau_y(ud, elem, node, 0.5)

        if hasattr(ud, 'rayleigh_forcing'):
                if ud.rayleigh_forcing:
                    ud.forcing_tcy, ud.forcing_tny = get_bottom_tau_y(ud, elem, node, 0.2, cutoff=0.1)

    n = 4.0           # eqn (12)
    A0 = 1.0e-3     # eqn (12)

    Msq = ud.Msq
    g = ud.gravity_strength[1]
    kappa = th.Gamma
    Omega = ud.coriolis_strength[2] / 2.0

    waveno = 1.0

    x = elem.x.reshape(-1,1)
    y = elem.y.reshape(1,-1)#[:,:ud.inbcy]

    Hrho = 1.0 / g
    kGam = 1.0 / Hrho * (1.0 / th.gamm - 1.0 / 2.0)       # eqn (8)
    A = A0 * bump(2.0 * y / (n * Hrho) - 1.0)             # eqn (12)

    use_hydrostate = False

    if use_hydrostate:
        # Use hydrostatically balanaced background
        hydrostatic_state(mpv, elem, node, th, ud)
        rhobar = mpv.HydroState.rho0.reshape(1,-1)
        Ybar = mpv.HydroState.Y0.reshape(1,-1)
        pibar = mpv.HydroState.p20.reshape(1,-1)
    else:
        # Use hydrostatic balance in Mark's notes
        Htheta = Hrho / kappa          
        rhobar = np.exp(-y / Hrho) 
        Ybar = np.exp(y / Htheta)  
        pibar = 1.0 / Ybar
        mpv.HydroState.rho0[...] = rhobar
        mpv.HydroState.Y0[...] = Ybar
        mpv.HydroState.S0[...] = 1.0 / Ybar #1.0 / ud.stratification(y)

        # pibar = (rhobar * Ybar * ud.Rg)**(th.gm1)
        # pibar = (rhobar * Ybar)**(th.gm1)

        # mpv.HydroState.rho0[:ud.inbcy] = rhobar
        # mpv.HydroState.Y0[:ud.inbcy] = Ybar

    # dimensionless Brunt-Väisälä frequency
    N = ud.t_ref * np.sqrt(ud.Nsq_ref) 
    # dimensionless speed of sound
    Cs = np.sqrt(th.gamm / Msq)  

    ud.u_wind_speed = 0.0 #-Cs

    # 1.0 / c_p in the notes
    # fac = 1.0 / (ud.cp_gas / (ud.h_ref**2 / ud.t_ref**2 / ud.T_ref))
    Gamma = 1.0 / Hrho * (1.0 / th.gamm - 0.5)

    # time shift of the initial solution
    ts = -0.5 / N * np.pi #/ ud.t_ref
    # ts = 50.0
    # ts = Cs * 200.0 / N 
    ts = 0.0

    # set up perturbation quantities
    # exp(-kGam * y)  / sqrt(rhobar) / Ybar = 1.0
    up = A0 * Ybar * np.cos(N / Cs * (waveno * x + Cs * ts))
    vp = 0.0
    wp = 0.0
    Yp = 0.0
    # th.Gamma * ud.Msq == 1 / dimensionless(c_p)
    fac = th.Gamma
    pi = A0 * Cs * fac * np.cos(N / Cs * (waveno * x + Cs * ts))

    # up = A / rhobar**0.5 * np.exp(-Gamma * y) * np.cos(N / Cs * (waveno * x + Cs * ts))
    # pi = A * Cs * fac / rhobar**0.5 / Ybar * np.exp(-Gamma * y) * np.cos(N / Cs * (waveno * x + Cs * ts))

    # up = up[:,::-1]
    # pi = pi[:,::-1]

    u = ud.u_wind_speed + up
    v = ud.v_wind_speed + vp
    w = ud.w_wind_speed + wp
    Y = Ybar + Yp

    # eqn (2.3)
    rho = (((pibar + Msq * pi))**th.gm1inv) / Y# / ud.Rg
    # rho = rhobar
    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhow[...] = rho * w
    Sol.rhoY[...] = rho * Y
    Sol.rhoX[...] = 0.0
    mpv.p2_cells[...] = pi


    # Sol.rho[:,:ud.inbcy] = rho
    # Sol.rhou[:,:ud.inbcy] = rho * u
    # Sol.rhov[:,:ud.inbcy] = rho * v
    # Sol.rhow[:,:ud.inbcy] = rho * w
    # Sol.rhoe[:,:ud.inbcy] = 0.0
    # Sol.rhoY[:,:ud.inbcy] = rho * Y
    # Sol.rhoX[:,:ud.inbcy] = 0.0
    # mpv.p2_cells[:,:ud.inbcy] = pi

    ###################################################
    # initialise nodal pi
    xn = node.x.reshape(-1,1)
    yn = node.y.reshape(1,-1)#[:,:ud.inbcy+1]

    # initialise nodal pressure
    Hrho_n = Hrho

    if use_hydrostate:
        # Use hydrostatically balanced background
        Ybar_n = mpv.HydroState_n.Y0.reshape(1,-1)
        rhobar_n = mpv.HydroState_n.rho0.reshape(1,-1)
    else:
        # Use hydrostatic balance from notes
        Ybar_n = np.exp(yn / Htheta)
        rhobar_n = np.exp(-yn / Hrho_n)
        mpv.HydroState_n.Y0[...] = Ybar_n
        mpv.HydroState_n.S0[...] = 1.0 / Ybar_n # 1.0 / ud.stratification(yn)

    An = A0 #* bump(2.0 * yn / (n * Hrho_n) - 1.0)

    pi_n = An * Cs * fac * np.cos(N / Cs * (waveno * xn + Cs * ts))
    # pi_n = An * Cs * fac / rhobar_n**0.5 / Ybar_n * np.exp(-Gamma * yn) * np.cos(N / Cs * (waveno * xn + Cs * ts))

    mpv.p2_nodes[...] = pi_n#[:,::-1]
    # mpv.p2_nodes[:,:ud.inbcy+1] = pi_n

    # rhobar_top = np.exp(-elem.y.reshape(1,-1)[:,ud.inbcy:] / Hrho)
    # Ybar_top = np.exp(elem.y.reshape(1,-1)[:,ud.inbcy:] / Htheta)
    # pibar_top = 1.0 / Ybar

    # rho = ((pibar + Msq * pi)**th.gm1inv) / Y

    # mpv.HydroState.rho0[ud.inbcy:] = rhobar_top
    # mpv.HydroState.Y0[ud.inbcy:] = Ybar_top

    # Sol.rho[:,ud.inbcy:] = Sol.rho[:,ud.inbcy-1].reshape(-1,1)
    # Sol.rhou[:,ud.inbcy:] = Sol.rhou[:,ud.inbcy-1].reshape(-1,1)
    # Sol.rhov[:,ud.inbcy:] = Sol.rhov[:,ud.inbcy-1].reshape(-1,1)
    # Sol.rhow[:,ud.inbcy:] = Sol.rhow[:,ud.inbcy-1].reshape(-1,1)
    # Sol.rhoY[:,ud.inbcy:] = Sol.rho[:,ud.inbcy:]*Ybar_top
    # mpv.p2_nodes[:,ud.inbcy+1:] = mpv.p2_nodes[:,ud.inbcy].reshape(-1,1)

    # cnst = 0.0
    # Sol.rho[:,ud.inbcy:] = cnst
    # Sol.rhou[:,ud.inbcy:] = cnst
    # Sol.rhov[:,ud.inbcy:] = cnst
    # Sol.rhow[:,ud.inbcy:] = cnst
    # Sol.rhoY[:,ud.inbcy:] = cnst
    # mpv.p2_nodes[:,ud.inbcy+1:] = cnst

    if ud.bdry_type[1] == 'RAYLEIGH':
        rayleigh_damping(Sol, mpv, ud, ud.tcy, elem, th)

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    return Sol
