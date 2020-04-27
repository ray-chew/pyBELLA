import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
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
        self.is_compressible = 0 # PSINC case
        self.is_ArakawaKonor = 0

        self.compressibility = 1.0
        self.acoustic_timestep = 0
        self.acoustic_order = 0
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)
        self.Msq = 0

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

        # self.xmin = - 0.5
        # self.xmax =   0.5
        # self.ymin = - 0.5
        # self.ymax =   0.5
        # self.zmin = - 0.5
        # self.zmax =   0.5

        self.xmin = 0.0
        self.xmax = 1.0
        self.ymin = 0.0
        self.ymax = 1.0
        self.zmin = 0.0
        self.zmax = 1.0

        self.wind_speed = 0.0
        self.wind_shear = -0.0
        self.hill_shape = HillShapes.AGNESI
        self.hill_height = 0.0
        self.hill_length_scale = 99999.9

        self.bdry_type_min = np.empty((3), dtype=object)
        self.bdry_type_max = np.empty((3), dtype=object)

        self.bdry_type_min[0] = BdryType.PERIODIC
        self.bdry_type_min[1] = BdryType.WALL
        self.bdry_type_min[2] = BdryType.PERIODIC
        self.bdry_type_max[0] = BdryType.PERIODIC
        self.bdry_type_max[1] = BdryType.WALL
        self.bdry_type_max[2] = BdryType.PERIODIC

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.WALL
        self.bdry_type[2] = BdryType.PERIODIC

        self.absorber = 0 # 0 == WRONG == FALSE 
        self.bottom_theta_bc = BottomBC.BOTTOM_BC_DEFAULT
        ##########################################
        # NUMERICS
        ##########################################

        self.time_integrator = TimeIntegrator.SI_MIDPT
        self.advec_time_integrator = TimeIntegrator.STRANG
        # self.CFL  = 0.9/2.0
        self.CFL = 0.95
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

        self.initial_projection = False
        self.initial_impl_Euler = False

        self.column_preconditionr = False
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 0.0

        self.tout = np.arange(0.0,1.001,0.005)

        self.stepmax = 10
        # self.stepmax = 20000

        self.output_base_name = "_travelling_vortex"
        self.output_name_psinc = "_low_mach_gravity_psinc"
        self.output_name_comp = "_low_mach_gravity_comp"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1])
        
        aux = '3D'
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
    igs = elem.igs
    igy = igs[1]

    igxn = node.igx
    igyn = node.igy

    i2, iy = [slice(None,)] * elem.ndim, [slice(None,)] * elem.ndim
    iy[1] = slice(igy,-igy)
    for dim in range(elem.ndim):
        i2[dim] = slice(igs[dim],-igs[dim])
    hi2 = np.copy(i2)
    hi2[1] = slice(None,)
    i2, iy, hi2 = tuple(i2), tuple(iy), tuple(hi2)

    u0 = 0.0 #* ud.wind_speed
    v0 = 0.0 #* ud.wind_speed
    w0 = 0.0 #* ud.wind_speed

    rotdir = 1.0

    p0 = 1.0
    a_rho = 1.0
    rho0 = a_rho * 0.5
    del_rho = a_rho * 0.5
    R0 = 0.4
    fac = 1. * 1024.0
    xc = 0.5
    yc = 0.5
    zc = 0.5
    aa = 1.0

    f = ud.coriolis_strength[0]

    if seed != None:
        np.random.seed(seed)
        xc = (np.random.random() - 0.5)/10.
        yc = (np.random.random() - 0.5)/10.
        zc = (np.random.random() - 0.5)/10.

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

    ccoe = np.zeros((19))
    ccoe[0] = 1.0 / 7.0
    ccoe[1] = -3.0 / 4.0
    ccoe[2] = 4.0 / 3.0
    ccoe[3] = -1.0 / 5.0
    ccoe[4] = -45.0 / 22.0
    ccoe[5] = 3.0 / 4.0
    ccoe[6] = 9.0 / 2.0
    ccoe[7] = -36.0 / 7.0
    ccoe[8] = -11.0 / 5.0
    ccoe[9] = 55.0 / 8.0
    ccoe[10] = -33.0 / 17.0
    ccoe[11] = -4.0
    ccoe[12] = 58.0 / 19.0
    ccoe[13] = 3.0 / 5.0
    ccoe[14] = -10.0 / 7.0
    ccoe[15] = 4.0 / 11.0
    ccoe[16] = 9.0 / 46.0
    ccoe[17] = -1.0 / 8.0
    ccoe[18] = 1.0 / 50.0

    xs = elem.x.reshape(-1,1,1)
    ys = elem.y[igy:-igy].reshape(1,-1,1)
    zs = elem.z.reshape(1,1,-1)
    xccs = np.zeros_like(xs)
    yccs = np.zeros_like(ys)
    zccs = np.zeros_like(zs)

    xcm = xc - (ud.xmax - ud.xmin)
    ycm = yc - (ud.ymax - ud.ymin)
    zcm = zc - (ud.zmax - ud.zmin)

    xccs[...] = xc * (np.abs(xs - xc) < np.abs(xs - xcm))
    xccs[...] += xcm * (np.abs(xs - xc) > np.abs(xs - xcm))

    yccs[...] = yc * (np.abs(ys - yc) < np.abs(ys - ycm))
    yccs[...] += ycm * (np.abs(ys - yc) > np.abs(ys - ycm))

    zccs[...] = zc * (np.abs(zs - zc) < np.abs(zs - zcm))
    zccs[...] += zcm * (np.abs(zs - zc) > np.abs(zs - zcm))

    # dxc = xc * (xs - xc)
    # dzc = zs - zc

    # rs = [(dxc)**2,(dzc)**2]
    rs = [(xs-xccs)**2,(zs-zccs)**2]

    rs = np.array(rs[:elem.ndim])
    r = np.sqrt(rs.sum())
    r = np.repeat(r,elem.icy,axis=1)[iy]
    
    uth = (rotdir * fac * (1.0 - r/R0)**6 * (r/R0)**6) * (r < R0)

    u = u0 + uth * (-(zs-zccs)/r)
    v = v0
    w = w0 + uth * (+(xs-xccs)/r)
    # p_hydro = mpv.HydroState.p0[igy:-igy]
    # rhoY = mpv.HydroState.rhoY0[igy:-igy]
    p_hydro = 1.0
    rhoY = np.ones_like(Sol.rhoY)

    rho = np.zeros_like(r)
    rho[...] += (rho0 + del_rho * (1. - (r/R0)**2)**6) * (r < R0)
    rho[...] += rho0 * (r > R0)

    rho, u, w = rho.squeeze(), u.squeeze(), w.squeeze()

    Sol.rho[iy] = rho
    Sol.rhou[iy] = rho * u
    Sol.rhov[iy] = rho * v
    Sol.rhow[iy] = rho * w

    dp2c = np.zeros_like((r))
    for ip in range(25):
        dp2c += fac * (a_rho * coe[ip] * ((r/R0)**(12+ip) - 1.0) * rotdir**2) * (r/R0 < 1.0)
    for ip in range(19):
        dp2c += f * ccoe[ip] * ((r/R0)**(7+ip) - 1.0) * (r/R0 < 1.0)

    p2c = np.copy(dp2c).squeeze()

    p = p0 + ud.Msq * fac * dp2c
    p = p.squeeze()
    # Sol.rhoY[iy] = p**th.gamminv
    Sol.rhoY = rhoY
    Sol.rhoe[iy] = ud.rhoe(rho,u,v,w,p,ud,th)

    # rhoY0c = mpv.HydroState.rhoY0[igy:-igy]

    # for dim in range(0,elem.ndim,2):
    #     rhoY0c = np.expand_dims(rhoY0c, dim)
    #     rhoY0c = np.repeat(rhoY0c, elem.sc[dim], axis=dim)

    # mpv.p2_cells[iy] = th.Gamma * fac**2 * np.divide(p2c, rhoY0c)

    # set_ghostcells_p2(mpv.p2_cells, elem, ud)

    xs = node.x[igs[0]:-igs[0]].reshape(-1,1,1)
    ys = node.y[igs[1]:-igs[1]].reshape(1,-1,1)
    zs = node.z[igs[2]:-igs[2]].reshape(1,1,-1)

    xccs = np.zeros_like(xs)
    yccs = np.zeros_like(ys)
    zccs = np.zeros_like(zs)

    yccs[np.where(np.abs(ys - yc) < np.abs(ys - ycm))] = yc
    yccs[np.where(np.abs(ys - yc) >= np.abs(ys - ycm))] = ycm

    xccs[np.where(np.abs(xs - xc) < np.abs(xs - xcm))] = xc
    xccs[np.where(np.abs(xs - xc) >= np.abs(xs - xcm))] = xcm

    zccs[np.where(np.abs(zs - zc) < np.abs(zs - zcm))] = zc
    zccs[np.where(np.abs(zs - zc) >= np.abs(zs - zcm))] = zcm

    # dxcn = xs - xc
    # dzcn = zs - zc

    # rs = [(dxcn)**2,(dzcn)**2]

    rs = [(xs-xccs)**2, (zs-zccs)**2]

    rs = np.array(rs[:elem.ndim])
    r = np.sqrt(rs.sum())
    r = np.repeat(r,node.icy,axis=1)[iy]
    
    for ip in range(25):
        mpv.p2_nodes[i2] += fac* (a_rho * coe[ip] * ((r/R0)**(12+ip) - 1.0) * rotdir**2)

    for ip in range(19):
        mpv.p2_nodes[i2] += f * ccoe[ip] * ((r/R0)**(7+ip) - 1.0)
    mpv.p2_nodes[i2] *= r/R0 < 1.0

    # rhoY0n = mpv.HydroState.rhoY0[igyn:-igyn+1]

    # for dim in range(0,node.ndim,2):
    #     rhoY0n = np.expand_dims(rhoY0n, dim)
    #     rhoY0n = np.repeat(rhoY0n, node.sc[dim], axis=dim)

    # mpv.p2_nodes[i2] = th.Gamma * fac**2 * np.divide(mpv.p2_nodes[i2] , rhoY0n[hi2])

    mpv.p2_nodes[i2] = th.Gamma * fac * mpv.p2_nodes[i2]

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)
    # set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    # Sol.rhoY[...] = 1.0
    # mpv.p2_nodes[...] = 0.0

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