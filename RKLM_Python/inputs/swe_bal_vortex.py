import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
from physics.hydrostatics import hydrostatic_state
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2
from physics.low_mach.second_projection import euler_backward_non_advective_impl_part

from scipy import signal
from scipy.interpolate import griddata

class UserData(object):
    NSPEC = 1

    grav = 0.0
    omega = 1.0#2.0 * np.pi #/ 100.0

    R_gas = 1.0
    R_vap = 1.0
    Q_vap = 1.0
    gamma = 2.0
    g0 = 9.81
    # g0 = 2.0

    viscm = 0.0
    viscbm = 0.0
    visct = 0.0
    viscbt = 0.0
    cond = 0.0

    h_ref = 1000.0
    d_ref = 10.0
    t_ref = 1000.0
    T_ref = 1.0 # can get rid.
    u_ref = h_ref / t_ref
    p_ref = 1.0
    rho_ref = p_ref / (R_gas * T_ref)

    Nsq_ref = 0.0

    i_gravity = np.zeros((3))
    i_coriolis = np.zeros((3))

    tout = np.zeros((2))

    def __init__(self):
        self.h_ref = self.h_ref
        self.d_ref = self.d_ref
        self.t_ref = self.t_ref
        self.T_ref = self.T_ref
        self.rho_ref = self.rho_ref
        self.u_ref = self.u_ref
        self.Nsq_ref = self.Nsq_ref
        self.g_ref = self.grav
        self.gamm = self.gamma
        self.Rg_over_Rv = self.R_gas / self.R_vap
        self.Q = self.Q_vap / (self.R_gas * self.T_ref)
        self.g0 = self.g0

        self.nspec = self.NSPEC

        self.is_nonhydrostatic = 1
        self.is_compressible = 1
        self.is_ArakawaKonor = 0

        self.compressibility = 1.0
        self.acoustic_timestep = 0
        self.acoustic_order = 0
        
        self.Msq = 2.0 * self.u_ref * self.u_ref / (self.d_ref * self.g0)
        self.g0 = self.g0 / (self.d_ref / self.t_ref**2)

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

        self.fac = 1E0
        self.xmin = - 0.5 * self.fac
        self.xmax =   0.5 * self.fac
        self.ymin = - 0.5
        self.ymax =   0.5
        self.zmin = - 0.5 * self.fac
        self.zmax =   0.5 * self.fac

        self.u_wind_speed = 1.0
        self.w_wind_speed = 1.0
        self.wind_shear = -0.0
        self.hill_shape = HillShapes.AGNESI
        self.hill_height = 0.0
        self.hill_length_scale = 99999.9

        self.bdry_type_min = np.empty((3), dtype=object)
        self.bdry_type_max = np.empty((3), dtype=object)

        self.bdry_type_min[0] = BdryType.PERIODIC
        self.bdry_type_min[1] = BdryType.PERIODIC
        self.bdry_type_min[2] = BdryType.PERIODIC
        self.bdry_type_max[0] = BdryType.PERIODIC
        self.bdry_type_max[1] = BdryType.PERIODIC
        self.bdry_type_max[2] = BdryType.PERIODIC

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.PERIODIC
        self.bdry_type[2] = BdryType.PERIODIC

        self.absorber = 0 # 0 == WRONG == FALSE 
        self.bottom_theta_bc = BottomBC.BOTTOM_BC_DEFAULT
        ##########################################
        # NUMERICS
        ##########################################

        self.time_integrator = TimeIntegrator.SI_MIDPT
        self.advec_time_integrator = TimeIntegrator.STRANG
        # self.CFL  = 0.
        self.CFL = 0.45

        self.dtfixed0 = 1.0 #/ self.t_ref
        self.dtfixed = 1.0 #/ self.t_ref

        self.inx = 64+1
        self.iny = 1+1
        self.inz = 64+1

        self.recovery_order = RecoveryOrder.SECOND
        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.kp = 0.0
        self.kz = 0.0
        self.km = 0.0
        self.kY = 0.0
        self.kZ = 0.0

        self.tol = 1.e-8
        self.max_iterations = 6000

        self.perturb_type = 'pos_perturb'
        self.blending_mean = 'rhoY' # 1.0, rhoY
        self.blending_conv = 'swe' #theta, rho, None

        self.initial_blending = False

        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.blending_weight = 0./16

        self.initial_projection = True
        self.initial_impl_Euler = False

        self.column_preconditionr = False
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 0.0

        # self.tout = np.arange(0, 3.0 + 0.01, 0.01)[1:]
        self.tout = np.arange(0, 1.0 + 0.01, 0.01)[1:]
        # self.tout = [1.0]
        # self.tout = np.arange(0, 3.0 + 0.01, 0.01)[1:]
        # self.tout = [1.0]
        # self.tout = [1E6]

        # self.tout = times.copy()

        # self.stepmax = 10
        self.stepmax = 100001

        self.output_base_name = "_swe_vortex"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1])
        
        # aux = 'wdawloc_pp_rhou_rhow_tra_ip_0.25_nonorm'
        aux = 'pp_tra_ip_nonorm'
        # aux = 'psinc_noib'
        # aux = 'comp_imbal_ib_full-12'
        # aux ='comp_debug_noib_imbal'
        # aux = 'debug_tra'
        # aux = 'comp_1.0_pp_tra_truth_ip'

        self.aux = aux
        # self.output_suffix = "_%i_%i_%i_%.1f_%s" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1],aux)

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_function

    def stratification_function(self, y):
        if type(y) == float:
            return 1.0
        else:
            return np.ones((y.shape))

    def rhoe_function(self,rho,u,v,w,p,ud,th):
        Msq = ud.compressibility * ud.Msq
        gm1inv = th.gm1inv

        return p * gm1inv + 0.5 * Msq * rho * (u**2 + v**2 + w**2)

def sol_init(Sol, mpv, elem, node, th, ud, seed=None):
    scale = ud.fac

    u0 = ud.u_wind_speed * scale
    v0 = 0.0 #* ud.wind_speed
    w0 = ud.w_wind_speed * scale

    rotdir = 1.0

    R0 = 0.4 * scale
    fac = 1. * 1024.0
    xc = 0.0
    zc = 0.0

    g = ud.g0
    f = ud.coriolis_strength[0]
    # g = 9.81
    H0 = 1.0

    if seed != None:
        np.random.seed(seed)
        # H0 += (np.random.random() - 0.5) * 0.1
        # print(seed, H0)
        xc += (np.random.random() - 0.5) * scale / 5.0
        zc += (np.random.random() - 0.5) * scale / 5.0
        print(seed, xc, zc)

    # used to generate truth

    if 'truth' in ud.aux:
        np.random.seed(2233)
        # H0 += (np.random.random() - 0.5) * 0.1
        # print(seed, H0
        # ud.H0 = H0

        xc += (np.random.random() - 0.5) * scale / 5.0
        zc += (np.random.random() - 0.5) * scale / 5.0
        print(seed, xc, zc)
        ud.xc = xc
        ud.zc = zc

    xcm = xc - np.sign(xc) * (ud.xmax - ud.xmin)
    zcm = zc - np.sign(zc) * (ud.zmax - ud.zmin)

    # xcm, zcm = 1.0, 1.0

    igs = elem.igs
    igy = igs[1]

    igxn = node.igx
    igzn = node.igz

    i2 = (slice(igs[0],-igs[0]),slice(igs[1],-igs[1]),slice(igs[2],-igs[2]))

    hydrostatic_state(mpv, elem, node, th, ud)

    xs = elem.x.reshape(-1,1,1)
    zs = elem.z.reshape(1,1,-1)

    # xs = np.linspace(-0.5,0.5,64).reshape(-1,1,1)
    # zs = np.linspace(-0.5,0.5,64).reshape(1,1,-1)
    
    xccs = np.zeros_like(xs)
    zccs = np.zeros_like(zs)

    xccs[...] = xc * (np.abs(xs - xc) < np.abs(xs - xcm))
    xccs[...] += xcm * (np.abs(xs - xc) > np.abs(xs - xcm))

    zccs[...] = zc * (np.abs(zs - zc) < np.abs(zs - zcm))
    zccs[...] += zcm * (np.abs(zs - zc) > np.abs(zs - zcm))

    r = np.sqrt((xs-xccs)**2 + (zs-zccs)**2)
    uth = (rotdir * fac * (1.0 - r/R0)**6 * (r/R0)**6) * (r < R0)

    u = u0 + uth * (-(zs-zccs)/r)
    w = w0 + uth * (+(xs-xccs)/r)
    v = v0

    coe = np.zeros((13))
    coe[0] = +1.0/12
    coe[1] = -12.0/13
    coe[2] = +33.0/7
    coe[3] = -44.0/3
    coe[4] = +495.0/16
    coe[5] = -792.0/17
    coe[6] = +154.0/3
    coe[7] = -792.0/19
    coe[8] = +99.0/4
    coe[9] = -220.0/21
    coe[10] = +3.0
    coe[11] = -12.0/23
    coe[12] = +1.0/24

    fcoe = np.zeros((7))
    fcoe[0] = +1.0/7
    fcoe[1] = -3.0/4
    fcoe[2] = +5.0/3
    fcoe[3] = -2.0
    fcoe[4] = +15.0/11
    fcoe[5] = -1.0/2
    fcoe[6] = +1.0/13

    rho = np.zeros_like(r)

    for i in range(12,24+1):
        rho[...] += fac**2 * coe[i-12] * (r/R0)**i * (r < R0)

    for i in range(7,13+1):
        rho[...] += f * fac * fcoe[i-7] * (r/R0)**i * (r < R0)

    p = np.copy(rho)
    rho *= ud.Msq / 2.0
    rho = (rho - rho.max()) * (r < R0)
    rho += H0

    if (ud.is_compressible):
        Sol.rho[:,igy:-igy,:] = rho
        Sol.rhou[:,igy:-igy,:] = rho * u
        Sol.rhov[:,igy:-igy,:] = rho * v
        Sol.rhow[:,igy:-igy,:] = rho * w
        Sol.rhoY[:,igy:-igy,:] = rho
    else:
        min_val = rho.mean()
        H0 = rho.mean()

        H1 = (rho - min_val)

        Sol.rho[:,igy:-igy,:] = H0
        Sol.rhou[:,igy:-igy,:] = H0 * u
        Sol.rhov[:,igy:-igy,:] = H0 * v
        Sol.rhow[:,igy:-igy,:] = H0 * w
        Sol.rhoY[:,igy:-igy,:] = H1

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    kernel = np.ones((2,2))
    kernel /= kernel.sum()
    if (ud.is_compressible):
        pn = (Sol.rhoY[:,igy,:] - H0) / ud.Msq
        pn = signal.convolve(pn, kernel, mode='valid')
    else:
        rhoY = Sol.rhoY[:,igy,:] / ud.Msq
        pn = signal.convolve(rhoY, kernel, mode='valid')
        Sol.rhoY[:,igy:-igy,:] = H0

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)
    
    pn = np.expand_dims(pn, axis=1)
    mpv.p2_nodes[1:-1,:,1:-1] = np.repeat(pn[...], node.icy, axis=1)
    mpv.p2_nodes[2:-2,2:-2,2:-2] -= mpv.p2_nodes[2:-2,2:-2,2:-2].mean(axis=(0,2),keepdims=True)

    # mpv.p2_nodes[...] = 1.0

    set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    if ud.initial_projection == True:
        is_compressible = np.copy(ud.is_compressible)
        compressibility = np.copy(ud.compressibility)
        ud.is_compressible = 0
        ud.compressibility = 0.0

        p2aux = np.copy(mpv.p2_nodes)

        Sol.rhou -= u0 * Sol.rho
        Sol.rhov -= v0 * Sol.rho
        Sol.rhow -= w0 * Sol.rho

        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, 0.0, ud.dtfixed, 0.5)
        # euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, 0.0, 2.5e-3, 0.5)

        mpv.p2_nodes[...] = p2aux
        mpv.dp2_nodes[...] = 0.0

        Sol.rhou += u0 * Sol.rho
        Sol.rhov += v0 * Sol.rho
        Sol.rhow += w0 * Sol.rho

        ud.is_compressible = is_compressible
        ud.compressibility = compressibility

    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)
