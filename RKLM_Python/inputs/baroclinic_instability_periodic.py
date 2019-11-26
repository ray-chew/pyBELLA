import numpy as np
import inputs.enum_bdry as enum_bdry
import management.enumerator as enumerator
import physics.hydrostatics as hydrostatics
import inputs.boundary as boundary
import management.variable as variable

from scipy import signal

class UserData(object):
    NSPEC = 1

    grav = 9.81
    omega = 1.0 * 0.0001
    R_gas = 287.4
    R_vap = 461.00
    Q_vap = 2.53e06
    gamma = 1.4

    h_ref = 1000
    t_ref = 100
    T_ref = 300.00
    p_ref = 1e5
    u_ref = h_ref / t_ref
    rho_ref = p_ref / (R_gas * T_ref)

    Nsq_ref = 0.00011772

    i_gravity = np.zeros((3))
    i_coriolis = np.zeros((3))

    def __init__(self):
        self.h_ref = self.h_ref
        self.t_ref = self.t_ref
        self.T_ref = self.T_ref
        self.p_ref = self.p_ref
        self.rho_ref = self.rho_ref
        self.u_ref = self.u_ref
        self.Nsq_ref = self.Nsq_ref
        self.g_ref = self.grav
        self.gamm = self.gamma
        self.Rg_over_Rv = self.R_gas / self.R_vap
        self.Q = self.Q_vap / (self.R_gas * self.T_ref)

        self.nspec = self.NSPEC

        self.is_ArakawaKonor = 0
        self.is_nonhydrostatic = 1
        self.is_compressible = 1
        self.acoustic_timestep = 0
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.gravity_strength = np.zeros((3))
        self.coriolis_strength = np.zeros((3))

        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.coriolis_strength[0] = self.omega * self.t_ref
        self.coriolis_strength[2] = self.omega * self.t_ref

        for i in range(3):
            if (self.gravity_strength[i] > np.finfo(np.float).eps) or (i == 1):
                self.i_gravity[i] = 1
                self.gravity_direction = 1

            if (self.coriolis_strength[i] > np.finfo(np.float).eps):
                self.i_coriolis[i] = 1

        self.xmin =   0.0e6 / self.h_ref
        self.xmax =  10.0e6 / self.h_ref
        self.ymin =   0.0e0 / self.h_ref
        self.ymax =  18.0e3 / self.h_ref
        self.zmin = - 8.0e6 / self.h_ref
        self.zmax =   8.0e6 / self.h_ref

        self.wind_speed = 0.0 * 20.0 / self.u_ref
        self.wind_shear = -0.0
        self.hill_shape = enumerator.HillShapes.AGNESI
        self.hill_height = 0.0 * 0.096447
        self.hill_length_scale = 0.1535

        self.bdry_type_min = np.empty((3), dtype=object)
        self.bdry_type_max = np.empty((3), dtype=object)

        self.bdry_type_min[0] = enum_bdry.BdryType.PERIODIC
        self.bdry_type_min[1] = enum_bdry.BdryType.WALL
        self.bdry_type_min[2] = enum_bdry.BdryType.PERIODIC
        self.bdry_type_max[0] = enum_bdry.BdryType.PERIODIC
        self.bdry_type_max[1] = enum_bdry.BdryType.WALL
        self.bdry_type_max[2] = enum_bdry.BdryType.PERIODIC

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = enum_bdry.BdryType.PERIODIC
        self.bdry_type[1] = enum_bdry.BdryType.WALL
        self.bdry_type[2] = enum_bdry.BdryType.PERIODIC

        self.absorber = 0 # 0 == WRONG == FALSE 
        self.bottom_theta_bc = enumerator.BottomBC.BOTTOM_BC_DEFAULT

        ##########################################
        # NUMERICS
        ##########################################
        self.time_integrator = enumerator.TimeIntegrator.SI_MIDPT
        self.advec_time_integrator = enumerator.TimeIntegrator.STRANG
        self.CFL  = 0.4
        self.dtfixed0 = 8.0
        self.dtfixed = 8.0

        # self.tips = TimeIntegratorParams()
        # SetTimeIntegratorParameters(self)

        self.inx = 65+1
        self.iny = 33+1
        self.inz = 65+1

        self.recovery_order = enumerator.RecoveryOrder.SECOND
        self.limiter_type_scalars = enumerator.LimiterType.NONE
        self.limiter_type_velocity = enumerator.LimiterType.NONE

        self.kp = 1.4
        self.kz = 1.4
        self.km = 1.4
        self.kY = 1.4
        self.kZ = 1.4

        self.continuous_blending = False
        self.no_of_pi_initial = 0
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        tol = 1.e-8
        self.flux_correction_precision = tol
        self.flux_correction_local_precision = tol
        self.second_projection_precision = tol
        self.second_projection_local_precision = tol
        self.flux_correction_max_iterations = 6000
        self.second_projection_max_iterations = 6000
        self.initial_projection = False

        self.column_preconditionr = 1
        # self.synchronize_nodal_pressure = False
        # self.synchronize_weight = 0.0

        self.eps_Machine = np.sqrt(np.finfo(np.float).eps)

        self.tout =  [2.5e1]
        # self.tout[1] = -1.0

        self.stepmax = 20000

        self.output_base_name = "_baroclinic_instability_periodic"
        self.output_name_psinc = "_low_mach_gravity_psinc"
        self.output_name_comp = "_low_mach_gravity_comp"
        self.output_suffix = ""

        self.ltrans = 525.0e3
        self.Temp0 = 273.0

        self.stratification = self.stratification_func
        self.zzero = self.zzero_func
        self.Thetae = self.Thetae_func
        self.Thetat = self.Thetat_func
        self.Thetas = self.Thetas_func
        self.Zeta = self.Zeta_func
        self.ht = self.ht_func

    def zzero_func(self,z):
        zmin = self.zmin
        zmax = self.zmax

        tmp = np.zeros_like(z)
        tmp[np.where(z >= 0.0)] = zmin + 3.0 * zmax
        tmp[np.where(z < 0.0)] = 3.0 * zmin + zmax
        return 0.25 * tmp
    
    def Thetae_func(self,y,z,inx):
        HTz = self.ht(z)

        # print(self.Thetat(y,z).shape)
        # print(self.Thetas(y,z).shape)
        # print((y < HTz))
        # print(~(y < HTz))

        tmp = (y < HTz) * self.Thetat(y,z) + ~(y < HTz) * self.Thetas(y,z)
        return np.repeat(tmp, inx, axis=0)
        # print(tmp.shape)
        # print(tmp[0])
        # print(tmp[1])

        # assert(0)

    def ht_func(self,z):
        dth = self.ltrans / self.h_ref      # [km]
        ht0 = 8.0e3 / self.h_ref            # [km]
        DTh = 30.0 / self.T_ref           # [K]
        Th0 = self.Temp0 / self.T_ref            # [K]
        Nsqs = 0.0245**2 * self.t_ref**2    # [s^{-2}]
        Nsqt = 0.01**2 * self.t_ref**2     # [s^{-2}]
        dNsq = Nsqs - Nsqt
        g = self.gravity_strength[1] / self.Msq

        dht = 0.5 * g * DTh / Th0 / dNsq
        fH0 = self.fH(z,dth)

        return (ht0 - dht * fH0)

    def fH(self, z, dth):
        z0 = self.zzero(z)
        return (np.sign(z) * self.F(z-z0,dth))

    @staticmethod
    def F(z, dth):
        eta = z / dth
        tmp = np.zeros_like(eta)
        tmp[...] = np.sin(np.pi * eta / 2.0)
        tmp[np.where(eta > 1.0)] = 1.0
        tmp[np.where(eta < -1.0)] = -1.0

        return tmp

    def Thetat_func(self,y,z):
        dth = self.ltrans / self.h_ref          # [km]
        DTh = 30.0 / self.T_ref                 # [K]
        Th0 = self.Temp0 / self.T_ref                # [K]

        Nsqt = 0.01**2 * self.t_ref**2          # [s^{-2}]
        g = self.gravity_strength[1] / self.Msq

        kappa = np.sign(z) * 70.0
        fHv = self.fH(z - kappa * y, dth)

        return Th0 * (1.0 + Nsqt * y / g) - 0.5 * DTh * fHv

    def Thetas_func(self,y,z):
        dth = self.ltrans / self.h_ref
        DTh = 30.0 / self.T_ref
        Th0 = self.Temp0 / self.T_ref

        Nsqs = 0.0245**2 * self.t_ref**2
        Nsqt = 0.01**2 * self.t_ref**2
        g = self.gravity_strength[1] / self.Msq

        # z0 = self.zzero(z)
        kappa = np.sign(z) * 70.0

        zeta = self.Zeta(y,z)

        HTz = self.ht(z)
        # H = 3.0e3 / self.h_ref
        # fHv = self.fH(z - kappa * HTz, dth)

        return Th0 * (1.0 + Nsqs * y / g) - Th0 * (Nsqs - Nsqt) * HTz * zeta / g - zeta * 0.5 * DTh * self.fH(z - kappa * HTz, dth) * (1.0 + 5.0 * (HTz - y) / HTz)


    def Zeta_func(self,y,z):
        H = 24.0e3 / self.h_ref
        return np.sin(0.5 * np.pi * (H - y) / (H - self.ht(z)))

    def stratification_func(self,y):
        Nsq = self.Nsq_ref * self.t_ref * self.t_ref
        g = self.gravity_strength[1] / self.Msq

        return np.exp(Nsq * y / g)




def sol_init(Sol, mpv, elem, node, th, ud, seed=None):
    xc = 0.5 * (ud.xmax + ud.xmin)
    yc = 9.0e3 / ud.h_ref
    wxz = 500.0e3 / ud.h_ref
    wy = 2.0e3 / ud.h_ref
    u0 = 0.0
    v0 = 0.0
    w0 = 0.0
    delth = 3.0 / ud.T_ref
    Ginv = th.Gammainv
    f = ud.coriolis_strength[0]
    dz = elem.dz
    dy = elem.dy

    icx = elem.icx
    icy = elem.icy
    icz = elem.icz

    igx = elem.igx
    igy = elem.igy
    igz = elem.igz

    icxn = node.icx
    icyn = node.icy

    g = ud.gravity_strength[1]

    hydrostatics.hydrostatic_state(mpv, elem, node, th, ud)

    HySt = variable.States([icxn,icyn], ud)
    HyStn = variable.States([icxn,icyn], ud)

    x = elem.x.reshape(-1,1,1)
    y = elem.y.reshape(1,-1,1)
    z = node.z.reshape(1,1,-1)

    zc = ud.zzero(z)
    r = np.sqrt( ((x-xc)**2 + (z-zc)**2) / (wxz**2) + (y-yc)**2 / (wy**2))
    r[np.where(r > 1.0)] = 1.0

    thp = delth * np.cos(0.5 * np.pi * r)**2.0
    the = ud.Thetae(y,z,elem.icx)

    Y = the + thp
    
    xn = node.x.reshape(-1,1,1)
    yn = node.y.reshape(1,-1,1)
    zn = node.z.reshape(1,1,-1)
    zc = ud.zzero(zn)
    r = np.sqrt( ((xn-xc)**2 + (zn-zc)**2) / (wxz**2) + (yn-yc)**2 / (wy**2) )
    r[np.where(r > 1.0)] = 1.0

    thpn = delth * np.cos(0.5 * np.pi * r)**2.0
    then = ud.Thetae(yn,zn,node.icx)
    Yn = then + thpn

    for k in range(icz):
        hydrostatics.hydrostatic_column(HySt, HyStn, Y[:,:,k], Yn[:-1,:-1,k], elem, node, th, ud)
        mpv.p2_nodes[:,:,k] = HyStn.p20[...]

    kernel = np.array([[[1.,1.],[1.,1.]], [[1.,1.],[1.,1.]]])
    kernel /= kernel.sum()

    p2 = signal.fftconvolve(mpv.p2_nodes, kernel, mode='valid') * ud.Msq

    # print(p2)
    # assert(0)
    P = p2**th.gm1inv
    p = p2**th.Gammainv

    kernel_dp2dy = np.array([[[1.,1.],[-1.,-1.]], [[1.,1.],[-1.,-1.]]])
    # kernel_dp2dy /= kernel_dp2dy.sum()

    dp2dy = 0.25 * signal.fftconvolve(mpv.p2_nodes, kernel_dp2dy, mode='valid') * ud.Msq / dy

    rhoY = P
    rho = -th.Gammainv * P * dp2dy / g

    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u0
    Sol.rhov[...] = rho * v0
    Sol.rhow[...] = rho * w0
    Sol.rhoY[...] = rhoY

    mpv.p2_cells[...] = p2 / ud.Msq
    Sol.rhoX[...] = rho / (P / rho)

    theta = Sol.rhoY / Sol.rho
    
    kernel_dpidz = np.array([[[-1.,-1.],[-1.,-1.]], [[1.,1.],[1.,1.]]])
    # kernel_dpidz /= kernel_dpidz.sum()

    dpidz = 0.25 * signal.fftconvolve(mpv.p2_nodes, kernel_dpidz, mode='valid') / dz
    u = -theta * Ginv * dpidz / f
    Sol.rhou = Sol.rho * u

    ud.nonhydrostasy = 0.0
    ud.compressibility = 0.0

    boundary.set_explicit_boundary_data(Sol, elem, ud, th, mpv)

    Ybar = Sol.rhoY / Sol.rho
    Ybar = np.mean(Ybar, axis=(0,2))

    mpv.HydroState.Y0[...] = Ybar
    mpv.HydroState.S0[...] = 1.0 / Ybar

    return Sol


    