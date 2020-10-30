import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
# from discretization.time_discretization import SetTimeIntegratorParameters
# from physics.gas_dynamics.explicit import TimeIntegratorParams
from physics.hydrostatics import hydrostatic_state
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2
from physics.low_mach.second_projection import euler_backward_non_advective_impl_part

class UserData(object):
    NSPEC = 1

    grav = 0.0
    omega = 0.0

    R_gas = 287.4
    R_vap = 461.0
    Q_vap = 2.53e+06
    gamma = 1.4

    viscm = 0.0
    viscbm = 0.0
    visct = 0.0
    viscbt = 0.0
    cond = 0.0

    h_ref = 10000
    t_ref = 100
    T_ref = 300.00
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

        self.xmin = - 0.5
        self.xmax =   0.5
        self.ymin = - 0.5
        self.ymax =   0.5
        self.zmin = - 0.5
        self.zmax =   0.5

        self.wind_speed = 0.0
        self.u_wind_speed = 1.0
        self.v_wind_speed = 1.0
        self.w_wind_speed = 0.0
        self.wind_shear = -0.0
        self.hill_shape = HillShapes.AGNESI
        self.hill_height = 0.0
        self.hill_length_scale = 99999.9

        # self.xmin = - 5000. / self.h_ref
        # self.xmax =   5000. / self.h_ref
        # self.ymin = - 5000. / self.h_ref / 8.0
        # self.ymax =   5000. / self.h_ref / 8.0
        # self.zmin = - 5000. / self.h_ref
        # self.zmax =   5000. / self.h_ref

        # self.wind_speed = 1.0 * 10.0 / self.u_ref
        # self.wind_shear = -0.0
        # self.hill_height = 0.0
        # self.hill_length_scale = 99999.9

        self.bdry_type_min = np.empty((3), dtype=object)
        self.bdry_type_max = np.empty((3), dtype=object)

        self.bdry_type_min[0] = BdryType.PERIODIC
        self.bdry_type_min[1] = BdryType.PERIODIC
        self.bdry_type_min[2] = BdryType.WALL
        self.bdry_type_max[0] = BdryType.PERIODIC
        self.bdry_type_max[1] = BdryType.PERIODIC
        self.bdry_type_max[2] = BdryType.WALL

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.PERIODIC
        self.bdry_type[2] = BdryType.WALL

        self.absorber = 0 # 0 == WRONG == FALSE 
        self.bottom_theta_bc = BottomBC.BOTTOM_BC_DEFAULT
        ##########################################
        # NUMERICS
        ##########################################

        self.time_integrator = TimeIntegrator.SI_MIDPT
        self.advec_time_integrator = TimeIntegrator.STRANG
        self.CFL  = 0.9/2.0
        # self.CFL = 0.95
        # self.dtfixed0 = 2.1 * 1.200930e-2
        # self.dtfixed = 2.1 * 1.200930e-2
        self.dtfixed = 0.01
        self.dtfixed0 = 0.01

        self.inx = 64+1
        self.iny = 64+1
        self.inz = 1

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
        self.blending_conv = 'rho' #theta, rho

        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.blending_weight = 10./16

        self.initial_blending = True

        self.initial_projection = False
        self.initial_impl_Euler = False

        self.column_preconditionr = False
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 0.0

        # self.tout = np.arange(0,11,1)/10.
        # self.tout = np.linspace(0,1.0,num=21)
        # self.tout = np.linspace(0.0,10.0,101)
        # self.tout = np.linspace(0.0,10.0,101*2-1)
        # self.tout = [0.5]
        # self.tout = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18,0.20]
        # self.tout = [0.1]
        # self.tout[0] =  1.0
        # self.tout[1] = -1.0
        self.tout = np.arange(0.0,3.01,0.01)[1:]
        # self.tout = np.arange(0.0,6.05,0.05)
        # self.tout = [1.0]


        # self.tout = times.copy()

        # self.stepmax = 10
        self.stepmax = 20000

        self.output_base_name = "_travelling_vortex"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.tout[-1])
        
        aux = 'ib_davg_noip'
        aux = 'ib_half-8_full-8'
        aux = 'ib-8_full-10'
        aux = 'comp_bal_noib'
        aux = 'comp_imbal_ib-16'
        # aux = 'psinc_noib'
        # aux = 'comp_1.0_pp_tra_truth'
        aux = 'wdawloc_pp_all_tra_0.25_nonorm'
        # aux = 'noda_pp'
        # aux = 'comp_debug_ib'
        self.aux = aux


        self.output_suffix = "_%i_%i_%.1f_%s" %(self.inx-1,self.iny-1,self.tout[-1],aux)

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
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = 0.0

    rotdir = 1.0

    p0 = 1.0
    a_rho = 1.0
    rho0 = a_rho * 0.5
    del_rho = a_rho * 0.5
    R0 = 0.4
    fac = 1. * 1024.0
    xc = 0.0
    yc = 0.0

    if seed != None and ud.perturb_type == 'pos_perturb':
        np.random.seed(seed)
        xc += (np.random.random() - 0.5) / 5.0
        yc += (np.random.random() - 0.5) / 5.0
        print(seed, xc, yc)

    if 'truth' in ud.aux:
        np.random.seed(2233)
        xc += (np.random.random() - 0.5) / 5.0
        yc += (np.random.random() - 0.5) / 5.0
        print(seed, xc, yc)
        ud.xc = xc
        ud.yc = yc

    xcm = xc - np.sign(xc) * (ud.xmax - ud.xmin)
    ycm = yc - np.sign(yc) * (ud.ymax - ud.ymin)

    igs = elem.igs
    igy = igs[1]

    igxn = node.igx
    igyn = node.igy

    hydrostatic_state(mpv, elem, node, th, ud)

    coe = np.zeros((25))
    coe[0]  =     1.0 / 12.0
    coe[1]  = -  12.0 / 13.0
    coe[2]  =     9.0 /  2.0
    coe[3]  = - 184.0 / 15.0
    coe[4]  =   609.0 / 32.0
    coe[5]  = - 222.0 / 17.0
    coe[6]  = -  38.0 /  9.0
    coe[7]  =    54.0 / 19.0
    coe[8]  =   783.0 / 20.0
    coe[9]  = - 558.0 /  7.0
    coe[10] =  1053.0 / 22.0
    coe[11] =  1014.0 / 23.0
    coe[12] = -1473.0 / 16.0
    coe[13] =   204.0 /  5.0
    coe[14] =   510.0 / 13.0
    coe[15] = -1564.0 / 27.0
    coe[16] =   153.0 /  8.0
    coe[17] =   450.0 / 29.0
    coe[18] = - 269.0 / 15.0
    coe[19] =   174.0 / 31.0
    coe[20] =    57.0 / 32.0
    coe[21] = -  74.0 / 33.0
    coe[22] =    15.0 / 17.0
    coe[23] = -   6.0 / 35.0
    coe[24] =     1.0 / 72.0

    xs = elem.x.reshape(-1,1)
    ys = elem.y[igy:-igy].reshape(1,-1)
    xccs = np.zeros_like(xs)
    yccs = np.zeros_like(ys)

    xccs[...] = xc * (np.abs(xs - xc) < np.abs(xs - xcm))
    xccs[...] += xcm * (np.abs(xs - xc) > np.abs(xs - xcm))

    yccs[...] = yc * (np.abs(ys - yc) < np.abs(ys - ycm))
    yccs[...] += ycm * (np.abs(ys - yc) > np.abs(ys - ycm))

    r = np.sqrt((xs-xccs)**2 + (ys-yccs)**2)
    uth = (rotdir * fac * (1.0 - r/R0)**6 * (r/R0)**6) * (r < R0)

    u = u0 + uth * (-(ys-yccs)/r)
    v = v0 + uth * (+(xs-xccs)/r)
    w = w0
    p_hydro = mpv.HydroState.p0[igy:-igy]
    rhoY = mpv.HydroState.rhoY0[igy:-igy]

    rho = np.zeros_like(r)
    rho[...] += (rho0 + del_rho * (1. - (r/R0)**2)**6) * (r < R0)
    rho[...] += rho0 * (r > R0)

    if seed != None and ud.perturb_type == 'fmp_perturb':
        np.random.seed(seed)
        rhop = np.fft.fft2(rho)
        rhop00 = np.copy(rhop[0,0])

        rhop = np.fft.fftshift(rhop)
        slc = (slice(int(elem.icx/2-1),int(elem.icx/2+2)),slice(int(elem.icy/2-1),int(elem.icy/2+2)))

        rhop[slc] += 30.0 * np.random.random(rhop[slc].shape)
        rhop = np.fft.ifftshift(rhop)
        rhop[0,0] = rhop00
        rho = np.fft.ifft2(rhop).real

    dp2c = np.zeros_like((r))
    for ip in range(25):
        dp2c += (a_rho * coe[ip] * ((r/R0)**(12+ip) - 1.0) * rotdir**2) * (r/R0 < 1.0)

    p2c = np.copy(dp2c)

    Sol.rho[:,igy:-igy] = rho
    Sol.rhou[:,igy:-igy] = rho * u
    Sol.rhov[:,igy:-igy] = rho * v
    Sol.rhow[:,igy:-igy] = rho * w

    if (ud.is_compressible) :
        p = p0 + ud.Msq * fac**2 * dp2c
        Sol.rhoY[:,igy:-igy] = p**th.gamminv
        Sol.rhoe[:,igy:-igy] = ud.rhoe(rho,u,v,w,p,ud,th)
    else:
        Sol.rhoe[:,igy:-igy] = ud.rhoe(rho,u,v,w,p_hydro,ud,th)
        Sol.rhoY[:,igy:-igy] = rhoY

    # mpv.p2_cells[:,igy:-igy] = th.Gamma * fac**2 * np.divide(p2c, mpv.HydroState.rhoY0[igy:-igy].T)

    mpv.p2_cells[:,igy:-igy] = th.Gamma * fac**2 * np.divide(p2c, mpv.HydroState.rhoY0[igy:-igy])    

    set_ghostcells_p2(mpv.p2_cells, elem, ud)

    xs = node.x[igxn:-igxn].reshape(-1,1)
    ys = node.y[igyn:-igyn].reshape(1,-1)
    # z = node.z
    xccs = np.zeros_like(xs)
    yccs = np.zeros_like(ys)

    yccs[np.where(np.abs(ys - yc) < np.abs(ys - ycm))] = yc
    yccs[np.where(np.abs(ys - yc) > np.abs(ys - ycm))] = ycm

    xccs[np.where(np.abs(xs - xc) < np.abs(xs - xcm))] = xc
    xccs[np.where(np.abs(xs - xc) > np.abs(xs - xcm))] = xcm

    r = np.sqrt((xs-xccs)**2 + (ys-yccs)**2)
    
    for ip in range(25):
        mpv.p2_nodes[igxn:-igxn,igyn:-igyn] += coe[ip] * ((r/R0)**(12+ip) - 1.0) * rotdir**2
    mpv.p2_nodes[igxn:-igxn,igyn:-igyn] *= r/R0 < 1.0

    mpv.p2_nodes[igxn:-igxn,igyn:-igyn] = th.Gamma * fac**2 * np.divide(mpv.p2_nodes[igxn:-igxn,igyn:-igyn] , mpv.HydroState.rhoY0[igyn:-igyn+1])
    # mpv.p2_nodes -= mpv.p2_nodes.mean()
    # mpv.p2_nodes[igxn:-igxn,igyn:-igyn] = th.Gamma * fac**2 * np.divide(mpv.p2_nodes[igxn:-igxn,igyn:-igyn] , mpv.HydroState.rhoY0[0,igyn:-igyn+1])

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    # ud.is_nonhydrostatic = 1
    # ud.is_compressible = 0
    ud.compressibility = float(ud.is_compressible)

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)
    # set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    # Sol.rhoY[...] = 1.0
    # mpv.p2_nodes[...] = 0.0

    # from scipy import signal
    # p2n = mpv.p2_nodes - mpv.p2_nodes.mean()
    # p2n -= p2n.min()
    # rhoY = (p2n + 1.0)**th.gm1inv - 1.0
    # rhoY = signal.fftconvolve(rhoY,np.ones([2] * rhoY.ndim),mode='valid') * 0.25
    # Sol.rhoY[...] = rhoY

    if ud.initial_projection == True:
        is_compressible = np.copy(ud.is_compressible)
        compressibility = np.copy(ud.compressibility)
        ud.is_compressible = 0
        ud.compressibility = 0.0

        p2aux = np.copy(mpv.p2_nodes)

        Sol.rhou -= u0 * Sol.rho
        Sol.rhov -= v0 * Sol.rho

        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, 0.0, ud.dtfixed, 0.5)
        # euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, 0.0, 2.5e-3, 0.5)

        mpv.p2_nodes[...] = p2aux
        mpv.dp2_nodes[...] = 0.0

        Sol.rhou += u0 * Sol.rho
        Sol.rhov += v0 * Sol.rho

        ud.is_compressible = is_compressible
        ud.compressibility = compressibility

    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)