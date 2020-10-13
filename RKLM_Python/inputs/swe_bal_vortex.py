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
    omega = 0.0#2.0 * np.pi #/ 100.0

    R_gas = 1.0
    R_vap = 461.0
    Q_vap = 2.53e+06
    gamma = 2.0

    viscm = 0.0
    viscbm = 0.0
    visct = 0.0
    viscbt = 0.0
    cond = 0.0

    h_ref = 1.0
    t_ref = 1.0
    T_ref = 1.0
    p_ref = 1e+5
    u_ref = h_ref / t_ref
    rho_ref = p_ref / (R_gas * T_ref)

    Nsq_ref = 0.0

    i_gravity = np.zeros((3))
    i_coriolis = np.zeros((3))

    tout = np.zeros((2))

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
        self.coriolis_strength[0] = self.omega * self.t_ref
        self.coriolis_strength[2] = self.omega * self.t_ref

        for i in range(3):
            if (self.gravity_strength[i] > np.finfo(np.float).eps) or (i == 1):
                self.i_gravity[i] = 1
                self.gravity_direction = 1

            if (self.coriolis_strength[i] > np.finfo(np.float).eps):
                self.i_coriolis[i] = 1

        self.fac = 1E0
        self.xmin = - 0.5 * self.fac#/ self.h_ref
        self.xmax =   0.5 * self.fac#/ self.h_ref
        self.ymin = - 0.5
        self.ymax =   0.5
        self.zmin = - 0.5 * self.fac#/ self.h_ref
        self.zmax =   0.5 * self.fac#/ self.h_ref

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
        # self.CFL  = 0.9
        self.CFL = 0.45
        self.dtfixed0 = 0.005 #/ self.t_ref
        self.dtfixed = 0.005 #/ self.t_ref

        self.dtfixed0 = 1.0 / self.t_ref
        self.dtfixed = 1.0 / self.t_ref

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

        self.continuous_blending = True
        self.no_of_pi_initial = 2
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.blending_weight = 16./16

        self.initial_projection = False
        self.initial_impl_Euler = False

        self.column_preconditionr = False
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 0.0

        self.tout = np.arange(0, 1.0 + 0.01, 0.01)[1:]
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
        
        aux = 'wdawloc_pp_rhou_rhow_tra'
        # aux = 'debug_psinc_1E3'
        # aux = 'comp_test_0'
        # aux = 'noda_pp'
        # aux = 'bld_test'
        # aux = 'comp_1.0_corr_1.0'
        # aux = 'debug_vortparam_lake_tra_corr_2pi'
        # aux = 'debug_H0'
        # aux = 'comp_1.0_tra_truth'
        # aux = 'comp_1.0_pp_tra_truth'
        # aux = 'get_initial_perturb'

        self.aux = aux
        self.output_suffix = "_%i_%i_%i_%.1f_%s" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1],aux)

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

    u0 = ud.u_wind_speed * scale#0.0 #* ud.wind_speed
    v0 = 0.0 #* ud.wind_speed
    w0 = ud.w_wind_speed * scale#0.0 *

    rotdir = 1.0

    R0 = 0.4 * scale
    fac = 1. * 1024.0
    xc = 0.0
    zc = 0.0

    g = 9.81
    H0 = 1.0
    # ud.Msq = (th.gamm / g)**2.0
    # ud.Msq = 1.0/(g)

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
        # print(seed, H0)
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

    rho = np.zeros_like(r)
    Frsq = (uth.max()**2 + uth.max()**2) / (g * H0)
    Frsq = 1.0 / (g * 1.0)
    # Frsq = 1.0 / (g * H0)

    for i in range(12,24+1):
        rho[...] += fac**2 * coe[i-12] * (r/R0)**i * (r < R0)

    rho *= Frsq * scale #* (r < R0) * scale
    rho = (rho - rho.max()) * (r < R0) #* H0
    rho += H0 #* (r < R0)

    if (ud.is_compressible):
        p = g / 2.0 * rho**2

        Sol.rho[:,igy:-igy,:] = rho #* 1.2
        Sol.rhou[:,igy:-igy,:] = rho * u
        Sol.rhov[:,igy:-igy,:] = rho * v
        Sol.rhow[:,igy:-igy,:] = rho * w

        Sol.rhoY[:,igy:-igy,:] = p**th.gamminv

        # Sol.rho[i2] = rho
        # Sol.rhou[i2] = rho * u
        # Sol.rhov[i2] = rho * v
        # Sol.rhow[i2] = rho * w
        # Sol.rhoY[i2] = p**th.gamminv
    else:
        min_val = H0

        # rho_diff = rho.max() - rho.min()
        rho_diff = H0 - rho.min()

        H1 = (rho - min_val) / rho_diff

        p = g / 2.0 * H0**2
        H1 = g / 2.0 * (H1)**2

        Sol.rho[:,igy:-igy,:] = H0
        Sol.rhou[:,igy:-igy,:] = H0 * u
        Sol.rhov[:,igy:-igy,:] = H0 * v
        Sol.rhow[:,igy:-igy,:] = H0 * w
        Sol.rhoY[:,igy:-igy,:] = H1**th.gamminv

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    kernel = np.ones((2,2))
    kernel /= kernel.sum()
    if (ud.is_compressible):
        pn = signal.convolve(Sol.rhoY[:,igy,:], kernel, mode='valid')
    else:
        pn = signal.convolve((Sol.rhoY[:,igy,:]), kernel, mode='valid')
        Sol.rhoY[:,igy:-igy,:] = p**th.gamminv
        set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)
    pn = np.expand_dims(pn, 1)
    pn = np.repeat(pn, node.icy, axis=1)

    mpv.p2_nodes[1:-1,:,1:-1] = pn
    set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    pn = np.expand_dims(mpv.p2_nodes[:,igy,:], 1)
    mpv.p2_nodes[...] = np.repeat(pn[...], node.icy, axis=1)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)