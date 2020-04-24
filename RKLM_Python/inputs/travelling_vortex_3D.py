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
    omega = 1.0

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
        # self.Msq = 1e-16

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

        self.wind_speed = 1.0
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
        self.CFL  = 0.9/2.0
        # self.CFL = 0.95
        self.dtfixed0 = 2.1 * 1.200930e-2
        self.dtfixed = 2.1 * 1.200930e-2

        self.inx = 64+1
        self.iny = 2+1
        # self.iny = 1+1
        self.inz = 64+1

        self.recovery_order = RecoveryOrder.SECOND
        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.kp = 0.0
        self.kz = 0.0
        self.km = 0.0
        self.kY = 0.0
        self.kZ = 0.0

        self.tol = 1.e-6
        self.max_iterations = 6000

        self.continuous_blending = False
        self.no_of_pi_initial = 0
        self.no_of_pi_transition = 10
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.initial_projection = True
        self.initial_impl_Euler = False

        self.column_preconditionr = False
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 0.0

        # self.tout = np.arange(0.0,1.1,0.25)
        self.tout = [1.0]

        # self.stepmax = 3
        self.stepmax = 20000

        self.output_base_name = "_travelling_vortex"
        self.output_name_psinc = "_low_mach_gravity_psinc"
        self.output_name_comp = "_low_mach_gravity_comp"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.tout[-1])
        
        aux = 'truth'
        # aux = 'very_low_mach'
        # aux = 'advected'
        # aux = 'even' if self.no_of_pi_initial % 2 == 0 else 'odd'
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
    u0 = 0.0 * ud.wind_speed
    v0 = 0.0 * ud.wind_speed
    w0 = 0.0 * ud.wind_speed

    rotdir = 1.0

    p0 = 1.0
    a_rho = 1.0
    rho0 = a_rho * 0.5
    del_rho = a_rho * 0.5
    R0 = 0.4
    fac = 1. * 1024.0
    xc = 0.0
    yc = 0.0
    zc = 0.0
    aa = 1.0

    if seed != None:
        np.random.seed(seed)
        xc = (np.random.random() - 0.5)/10.
        yc = (np.random.random() - 0.5)/10.
        zc = (np.random.random() - 0.5)/10.

    xcm = xc - (ud.xmax - ud.xmin)
    ycm = yc - (ud.ymax - ud.ymin)
    zcm = zc - (ud.zmax - ud.zmin)

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

    xs = elem.x.reshape(-1,1,1)
    ys = elem.y[igy:-igy].reshape(1,-1,1)
    zs = elem.z.reshape(1,1,-1)
    xccs = np.zeros_like(xs)
    yccs = np.zeros_like(ys)
    zccs = np.zeros_like(zs)

    xccs[...] = xc * (np.abs(xs - xc) < np.abs(xs - xcm))
    xccs[...] += xcm * (np.abs(xs - xc) > np.abs(xs - xcm))

    yccs[...] = yc * (np.abs(ys - yc) < np.abs(ys - ycm))
    yccs[...] += ycm * (np.abs(ys - yc) > np.abs(ys - ycm))

    zccs[...] = zc * (np.abs(zs - zc) < np.abs(zs - zcm))
    zccs[...] += zcm * (np.abs(zs - zc) > np.abs(zs - zcm))

    rs = [(xs-xccs)**2,(ys-yccs)**2,(zs-zccs)**2]

    rs = np.array(rs[:elem.ndim])
    r = np.sqrt(rs.sum())
    
    uth = (rotdir * fac * (1.0 - r/R0)**6 * (r/R0)**6) * (r < R0)

    u = u0 + uth * (-(ys-yccs)/r)
    # u = u0 + uth * (-(zs-zccs)/r)
    v = v0 + uth * (+(xs-xccs)/r)
    # w = w0 + uth * (+(xc-xccs)/r)
    w = w0
    p_hydro = mpv.HydroState.p0[igy:-igy]
    rhoY = mpv.HydroState.rhoY0[igy:-igy]

    rho = np.zeros_like(r)
    rho[...] += (rho0 + del_rho * (1. - (r/R0)**2)**6) * (r < R0)
    rho[...] += rho0 * (r > R0)

    dp2c = np.zeros_like((r))
    for ip in range(25):
        dp2c += (a_rho * coe[ip] * ((r/R0)**(12+ip) - 1.0) * rotdir**2) * (r/R0 < 1.0)

    p2c = np.copy(dp2c).squeeze()

    i2, iy = [slice(None,)] * elem.ndim, [slice(None,)] * elem.ndim
    iy[1] = slice(igy,-igy)
    for dim in range(elem.ndim):
        i2[dim] = slice(igs[dim],-igs[dim])
    hi2 = np.copy(i2)
    hi2[1] = slice(None,)
    i2, iy, hi2 = tuple(i2), tuple(iy), tuple(hi2)

    rho, u, v = rho.squeeze(), u.squeeze(), v.squeeze()

    Sol.rho[iy] = rho
    Sol.rhou[iy] = rho * u
    Sol.rhov[iy] = rho * v
    Sol.rhow[iy] = rho * w

    if (ud.is_compressible) :
        p = p0 + ud.Msq * fac**2 * dp2c
        p = p.squeeze()
        Sol.rhoY[iy] = p**th.gamminv
        # Y = np.zeros_like(Sol.rhoY)
        # Y[:,igy:-igy] = Sol.rhoY[:,igy:-igy] / Sol.rho[:,igy:-igy]
        # Sol.rhoX[:,igy:-igy] = Sol.rho[:,igy:-igy] / Y[:,igy:-igy]
        Sol.rhoe[iy] = ud.rhoe(rho,u,v,w,p,ud,th)
    else:
        Sol.rhoe[iy] = ud.rhoe(rho,u,v,w,p_hydro,ud,th)
        Sol.rhoY[iy] = rhoY

    # mpv.p2_cells[:,igy:-igy] = th.Gamma * fac**2 * np.divide(p2c, mpv.HydroState.rhoY0[igy:-igy].T)
    rhoY0c = mpv.HydroState.rhoY0[igy:-igy]

    for dim in range(0,elem.ndim,2):
        rhoY0c = np.expand_dims(rhoY0c, dim)
        rhoY0c = np.repeat(rhoY0c, elem.sc[dim], axis=dim)

    mpv.p2_cells[iy] = th.Gamma * fac**2 * np.divide(p2c, rhoY0c)

    set_ghostcells_p2(mpv.p2_cells, elem, ud)

    xs = node.x[igs[0]:-igs[0]].reshape(-1,1,1)
    ys = node.y[igs[1]:-igs[1]].reshape(1,-1,1)
    zs = node.z[igs[2]:-igs[2]].reshape(1,1,-1)

    xccs = np.zeros_like(xs)
    yccs = np.zeros_like(ys)
    zccs = np.zeros_like(zs)

    yccs[np.where(np.abs(ys - yc) < np.abs(ys - ycm))] = yc
    yccs[np.where(np.abs(ys - yc) > np.abs(ys - ycm))] = ycm

    xccs[np.where(np.abs(xs - xc) < np.abs(xs - xcm))] = xc
    xccs[np.where(np.abs(xs - xc) > np.abs(xs - xcm))] = xcm

    zccs[np.where(np.abs(zs - zc) < np.abs(zs - zcm))] = zc
    zccs[np.where(np.abs(zs - zc) > np.abs(zs - zcm))] = zcm

    # r = np.sqrt((xs-xccs)**2 + (ys-yccs)**2)
    rs = [(xs-xccs)**2,(ys-yccs)**2,(zs-zccs)**2]

    rs = np.array(rs[:elem.ndim])
    r = np.sqrt(rs.sum().squeeze())
    
    for ip in range(25):
        mpv.p2_nodes[i2] += coe[ip] * ((r/R0)**(12+ip) - 1.0) * rotdir**2
    mpv.p2_nodes[i2] *= r/R0 < 1.0

    rhoY0n = mpv.HydroState.rhoY0[igyn:-igyn+1]

    for dim in range(0,node.ndim,2):
        rhoY0n = np.expand_dims(rhoY0n, dim)
        rhoY0n = np.repeat(rhoY0n, node.sc[dim], axis=dim)

    mpv.p2_nodes[i2] = th.Gamma * fac**2 * np.divide(mpv.p2_nodes[i2] , rhoY0n[hi2])
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

        mpv.p2_nodes[...] = p2aux
        mpv.dp2_nodes[...] = 0.0

        Sol.rhou += u0 * Sol.rho
        Sol.rhov += v0 * Sol.rho

        ud.is_compressible = is_compressible
        ud.compressibility = compressibility

    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)