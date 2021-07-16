import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
from physics.hydrostatics import hydrostatic_state, hydrostatic_column, hydrostatic_initial_pressure
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2
from management.variable import States

class UserData(object):
    NSPEC = 1
    BOUY = 0

    grav = 9.80665
    omega = 1.0*0.000103126

    R_gas = 287.05
    R_vap = 461.0
    Q_vap = 2.53e+06
    gamma = 1.4

    T_ref = 250.00                              # [K]
    u_ref = 10.0                                # [m/s]
    p_ref = 1e+5                                # [Pa]
    h_ref = R_gas * T_ref / grav                # [m]
    t_ref = h_ref / u_ref                       # [s]
    rho_ref = p_ref / (R_gas * T_ref)

    Nsq_ref = ((gamma-1.0)/gamma) * grav * grav / (R_gas*T_ref)

    # planetary -> 160.0;  long-wave -> 20.0;  standard -> 1.0;
    scale_factor = 20.0

    i_gravity = np.zeros((3))
    i_coriolis = np.zeros((3))

    tout = np.zeros((2))

    def __init__(self):
        self.scale_factor = self.scale_factor

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

        self.is_nonhydrostatic = 0
        self.is_compressible = 1
        self.is_ArakawaKonor = 0

        self.compressibility = 0.0
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
                self.gravity_direction = i

            if (self.coriolis_strength[i] > np.finfo(np.float).eps):
                self.i_coriolis[i] = 1

        self.xmin = -0.5 * self.scale_factor * 300000.0 / self.h_ref
        self.xmax =  0.5 * self.scale_factor * 300000.0 / self.h_ref
        self.ymin =  0.0
        self.ymax =  10000.0 / self.h_ref
        self.zmin = -1.0
        self.zmax =  1.0

        self.u_wind_speed = 0.0 * 20.0 / self.u_ref
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.WALL
        self.bdry_type[2] = BdryType.WALL

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL  = 0.9
        self.dtfixed0 = 0.5 / self.t_ref * 50.0 * self.scale_factor
        self.dtfixed = 0.5 / self.t_ref * 50.0 * self.scale_factor
        # self.dtfixed0 = 0.5 / self.t_ref
        # self.dtfixed = 0.5 / self.t_ref

        self.inx = 301+1
        # self.inx = 1201+1
        self.iny = 20+1
        # self.iny = 80+1
        self.inz = 1

        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.initial_projection = False

        self.tout = np.array([8 * 180.0 * self.scale_factor / self.t_ref]) # 8 hrs
        # self.tout = np.append(self.tout, 0.17081780*4.0)
        # self.tout = np.sort(self.tout)
        # self.tout = np.arange(0,self.tout[0],900 / self.t_ref)[1:]

        self.tol = 1.e-8
        self.stepmax = 1000000
        self.max_iterations = 6000

        self.initial_blending = False
        self.no_of_pi_initial = 0
        self.no_of_hy_initial = 1

        self.continuous_blending = False
        self.blending_weight = 0./16
        self.blending_mean = 'rhoY' # 1.0, rhoY
        self.blending_conv = 'rho' #theta, rho
        self.blending_type = 'half' # half, full

        self.output_base_name = "_baldauf_brdar"
        if self.scale_factor < 10.0:
            self.scale_tag = "standard"
        elif self.scale_factor == 20.0:
            self.scale_tag = "long"
        elif self.scale_factor == 160.0:
            self.scale_tag = "planetary"
        else:
            assert(0, "scale factor unsupported")

        if self.is_nonhydrostatic and self.is_compressible:
            self.h_tag = 'nonhydro'
        elif self.is_nonhydrostatic and not self.is_compressible:
            self.h_tag = 'psinc'
        elif not self.is_nonhydrostatic:
            self.h_tag = 'hydro'

        aux = ''
        if len(aux) > 0:
            aux = self.scale_tag + "_" + self.h_tag + "_" + aux
        else:
            aux = self.scale_tag + "_" + self.h_tag
        self.aux = aux

        self.output_timesteps = False

        self.stratification = self.stratification_function
        self.molly = self.molly_function
        self.rhoe = self.rhoe_method
        
    def stratification_function(self, y):
        Nsq = self.Nsq_ref * self.t_ref * self.t_ref
        g = self.gravity_strength[1] / self.Msq

        return np.exp(Nsq * y / g)

    def molly_function(self, x):
        del0 = 0.25
        L = self.xmax - self.xmin
        xi_l = np.minimum(1.0, (x - self.xmin) / (del0 * L))
        xi_r = np.minimum(1.0, (self.xmax - x) / (del0 * L))

        return 0.5 * np.minimum(1.0 - np.cos(np.pi * xi_l), 1.0 - np.cos(np.pi * xi_r))

    @staticmethod
    def rhoe_method(rho,u,v,w,p,ud,th):
        Msq = ud.compressibility * ud.Msq

        gm1inv = th.gm1inv
        return (p * gm1inv + 0.5 * Msq * rho * (u*u + v*v + w*w))

def sol_init(Sol, mpv, elem, node, th, ud, seeds=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed
    delT = 0.01 / ud.T_ref
    xc = -0.0 * ud.scale_factor * 50.0e+3 / ud.h_ref
    a = ud.scale_factor * 5.0e+3 / ud.h_ref
    H = ud.ymax - ud.ymin

    hydrostatic_state(mpv, elem, node, th, ud)
    
    x = elem.x.reshape(-1,1)
    y = elem.y.reshape(1,-1)

    Tb = delT * ud.molly(x) * np.sin(np.pi * y / H) * np.exp(-(x-xc)**2/(a**2))
    # rhoprime = -np.exp(-0.5 * delta * y) * Tb
    xn = node.x[:-1].reshape(-1,1)
    yn = node.y[:-1].reshape(1,-1)

    Tbn = delT * ud.molly(xn) * np.sin(np.pi * yn / H) * np.exp(-(xn-xc)**2/(a**2))

    HySt = States(node.sc,ud)
    HyStn = States(node.sc,ud)

    Y = ud.stratification(y) + Tb
    Yn = ud.stratification(yn) + Tbn

    hydrostatic_column(HySt, HyStn, Y, Yn, elem, node, th, ud)

    # p = mpv.HydroState.p0.reshape(1,-1)
    # rho = mpv.HydroState.rho0.reshape(1,-1) + rhoprime
    # rhoY = p**th.gamminv

    xc_idx = slice(0,-1)
    yc_idx = slice(0,-1)
    c_idx = (xc_idx, yc_idx)

    if ud.is_compressible:
        p = HySt.p0[c_idx]
        rhoY = HySt.rhoY0[c_idx]
    else:
        p = mpv.HydroState.p0
        rhoY = mpv.HydroState.rhoY0

    u, v, w = u0, v0, w0

    rho = rhoY / Y
    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhow[...] = rho * w
    Sol.rhoe[...] = ud.rhoe(rho, u ,v, w, p, ud, th)
    Sol.rhoY[...] = rhoY

    mpv.p2_cells[...] = HySt.p20[c_idx]

    Sol.rhoX[...] = Sol.rho * (Sol.rho / Sol.rhoY - mpv.HydroState.S0.reshape(1,-1))

    mpv.p2_nodes[:,:] = HyStn.p20

    hydrostatic_initial_pressure(Sol,mpv,elem,node,ud,th)

    if 'imbal' in ud.aux:
        mpv.p2_nodes[...] = 0.0

    mpv.HyStn = HyStn

    ud.nonhydrostasy = 1.0 if ud.is_nonhydrostatic == 1 else 0.0
    ud.compressibility = 1.0 if ud.is_compressible == 1 else 0.0

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)