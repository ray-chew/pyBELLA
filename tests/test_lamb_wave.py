import numpy as np
from dycore.utils.enum_bdry import BdryType
from dycore.utils.options import LimiterType
from dycore.physics.hydrostatics import hydrostatic_state
from dycore.utils.boundary import set_explicit_boundary_data, get_tau_y, get_bottom_tau_y

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
        self.R_gas = self.R_gas
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
        self.ymax =   2.0
        self.zmin = - 1.0
        self.zmax =   1.0

        self.u_wind_speed = 0.0 * self.u_ref
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.WALL
        self.bdry_type[2] = BdryType.WALL
        self.LAMB_BDRY = False

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.9

        self.inx = 151+1
        self.iny = 15+1
        self.inz = 1

        self.dtfixed0 = 100.0 / self.t_ref
        self.dtfixed = self.dtfixed0
        
        self.do_advection = False
        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.tol = 1.e-30
        self.max_iterations = 10000

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


        self.tout = [36.0]
        # self.tout = np.arange(0,361,1.0)
        # self.tout = np.append(self.tout, [720.0])
        self.stepmax = 31
        self.output_timesteps = True

        self.autogen_fn = False

        self.output_base_name = "_lamb_wave"
        self.output_type = 'test'
        self.aux = ''
        self.output_suffix = "_%i_%i" %(self.inx-1,self.iny-1)

        self.diag = True
        self.diag_current_run = 'test_lamb_wave'
        self.diag_plot_compare = True

        self.stratification = self.stratification_wrapper
        self.rayleigh_bc = self.rayleigh_bc_function
        self.init_forcing = self.forcing

        self.rayleigh_forcing = True
        self.rayleigh_forcing_type = 'func' # func or file
        self.rayleigh_forcing_fn = 'output_mark_wave_ensemble=1_601_240_bottom_forcing_S16.h5'
        self.rayleigh_forcing_path = './output_mark_wave'
        

    def stratification_wrapper(self, dy):
        return lambda y : self.stratification_function(y, dy)

    def stratification_function(self, y, dy):
        g = self.gravity_strength[1]

        Gamma = (self.gamm - 1.0) / self.gamm
        Hex = 1.0 / (Gamma * g)
        pi_m = np.exp(-(y - 0.5 * dy) / Hex)
        pi_p = np.exp(-(y + 0.5 * dy) / Hex)
        
        Theta = - (Gamma * g * dy) / (pi_p - pi_m)

        return Theta

        # Nsq = self.Nsq_ref * self.t_ref**2
        # g = self.gravity_strength[1] / self.Msq

        # return np.exp(Nsq * y / g)

    @staticmethod
    def rayleigh_bc_function(ud):
        if ud.bdry_type[1] == BdryType.RAYLEIGH or ud.rayleigh_forcing == True:
            ud.inbcy = ud.iny - 1
            ud.iny0 = np.copy(ud.iny)
            ud.iny = ud.iny0 + int(3*ud.inbcy)

            # tentative workaround
            ud.bcy = ud.ymax
            ud.ymax += 3.0 * ud.bcy


    class forcing(object):
        def __init__(self, k, mu, Cs, F, N, Gamma, ampl, g, rhobar, Ybar, rhobar_n, Ybar_n, X, Y, Xn, Yn):
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
            matrix = -np.array([[0, self.F, 0, 1j*self.Cs*self.k], 
                            [-self.F, 0, -self.N, self.Cs*(self.mu+self.Gamma)], 
                            [0, self.N, 0, 0], 
                            [1j*self.Cs*self.k, self.Cs*(self.mu-self.Gamma), 0, 0]])
        
            self.T_matrix = matrix

        def eigenfunction(self, t, s, grid='c'):
            if grid == 'c':
                x, z = self.X, self.Y
            elif grid == 'n':
                x, z, = self.Xn, self.Yn
            
            # Compute eigenvalues and eigenvectors
            eigval, eigvec = np.linalg.eig( self.T_matrix )

            # Find index of eigenvalue 
            # with greatest real part aka the instability growth rate
            ind = np.argmax( np.real( eigval ) )

            # construct solution according to eq. 2.27 and 2.19
            exponentials = np.exp( 1j * self.k * x + self.mu * z 
                                + ( eigval[ind] ) * (t) + 1j * s * t )
            chi_u  = self.ampl * np.real( eigvec[0,ind] * exponentials )
            chi_w  = self.ampl * np.real( eigvec[1,ind] * exponentials )
            chi_th = self.ampl * np.real( eigvec[2,ind] * exponentials )
            chi_pi = self.ampl * np.real( eigvec[3,ind] * exponentials )

            self.arrs = ( chi_u, chi_w, chi_th, chi_pi )

        def dehatter(self, th, grid='c'):
            if grid == 'n':
                Ybar = self.Ybar_n
                oorhobarsqrt = self.oorhobarsqrt_n
            elif grid == 'c':
                Ybar = self.Ybar
                oorhobarsqrt = self.oorhobarsqrt

            chi_u, chi_v, chi_Y, chi_pi = self.arrs

            up = oorhobarsqrt * chi_u
            vp = oorhobarsqrt * chi_v
            Yp = oorhobarsqrt * self.N / self.g * Ybar * chi_Y
            pi_p = oorhobarsqrt * self.Cs / Ybar / th.Gammainv * chi_pi
            
            return up.T, vp.T, Yp.T, pi_p.T


def sol_init(Sol, mpv, elem, node, th, ud, seeds=None):

    if hasattr(ud, 'rayleigh_bdry_switch'):
        if ud.rayleigh_bdry_switch:
            ud.bdry_type[1] = BdryType.RAYLEIGH

    if ud.bdry_type[1] == BdryType.RAYLEIGH:
        ud.tcy, ud.tny = get_tau_y(ud, elem, node, 0.5)

    if ud.rayleigh_forcing:
        # ud.tcy, ud.tny = get_tau_y(ud, elem, node, 0.005)
        
        ud.forcing_tcy, ud.forcing_tny = get_bottom_tau_y(ud, elem, node, 0.2, cutoff=0.3)

    A0 = 1.0e-1 / ud.u_ref
    Msq = ud.Msq
    g = ud.gravity_strength[1] * ud.Rg

    x = elem.x.reshape(-1,1)
    y = elem.y.reshape(1,-1)
    X, Y = np.meshgrid(x, y)

    xn = node.x.reshape(-1,1)
    yn = node.y.reshape(1,-1)

    dy = np.diff(node.y)[0]

    Xn, Yn = np.meshgrid(xn, yn)

    ##################################################
    # Following Rupert's fix, reinitialise all background quantities
    # as derived from one quantity.
    ud.stratification = ud.stratification(dy)

    # Use hydrostatically balanced background
    hydrostatic_state(mpv, elem, node, th, ud)
    rhobar = mpv.HydroState.rho0.reshape(1,-1)
    Ybar = mpv.HydroState.Y0.reshape(1,-1)
    pibar = mpv.HydroState.p20.reshape(1,-1) * ud.Msq

    rhobar_n = mpv.HydroState_n.rho0.reshape(1,-1)
    Ybar_n = mpv.HydroState_n.Y0.reshape(1,-1)

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
    # if F == 0.0:
    #     F += 1e-15

    G =  np.sqrt( 9. / 40. )
    Gamma = G * N / Cs
    k = N / Cs 

    ud.rf_bot = ud.init_forcing(k, -Gamma, Cs, F, N, Gamma, A0, g, rhobar, Ybar, rhobar_n, Ybar_n, X, Y, Xn, Yn)
    ud.rf_bot.get_T_matrix()

    ud.u_wind_speed = 0.0

    ud.rf_bot.eigenfunction(0, 1)
    up, vp, Yp, pi_p = ud.rf_bot.dehatter(th)

    u = ud.u_wind_speed + up
    v = ud.v_wind_speed + vp
    w = ud.w_wind_speed
    Y = Ybar + Yp

    # rho = (((pibar + Msq * pi_p))**th.gm1inv) / Y
    rho = rhobar

    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhow[...] = rho * w
    Sol.rhoY[...] = rho * Y
    Sol.rhoX[...] = 0.0
    mpv.p2_cells[...] = pi_p

    ###################################################
    # initialise nodal pi
    ud.rf_bot.eigenfunction(0, 1, grid='n')
    _, _, _, pi_n = ud.rf_bot.dehatter(th, grid='n')

    mpv.p2_nodes[...] = pi_n

    # if ud.bdry_type[1] == 'RAYLEIGH':
    #     rayleigh_damping(Sol, mpv, ud, ud.tcy, elem, th)

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    if hasattr(ud, 'mixed_run'):
        if ud.mixed_run:
            ud.coriolis_strength[2] = 2.0 * 7.292 * 1e-5 * ud.t_ref

    if hasattr(ud, 'trad_forcing'):
        if ud.trad_forcing:
            ud.rf_bot.F = 0.0
            ud.rf_bot.get_T_matrix()

    return Sol