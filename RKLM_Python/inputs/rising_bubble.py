import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
# from discretization.time_discretization import SetTimeIntegratorParameters
# from physics.gas_dynamics.explicit import TimeIntegratorParams
from physics.hydrostatics import hydrostatic_state, hydrostatic_initial_pressure
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2
from management.variable import States

from physics.low_mach.second_projection import euler_backward_non_advective_impl_part

class UserData(object):
    NSPEC = 1
    BOUY = 0

    grav = 10.0             # [m/s^2]
    omega = 0.0 * 2.0 * np.pi * np.sin(0.25 * np.pi) / (24.0 * 3600.0)

    R_gas = 287.4
    R_vap = 461.0
    Q_vap = 2.53e+06
    gamma = 1.4

    h_ref = 10000.0         # [m]
    t_ref = 100.0           # [s]
    T_ref = 300.00          # [K]
    p_ref = 8.61 * 1e4      # [N/m^2]
    u_ref = h_ref / t_ref
    rho_ref = p_ref / (R_gas * T_ref)

    Nsq_ref = grav * 1.3e-05

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
        self.bouy = self.BOUY

        self.mol_trans = MolecularTransport.NO_MOLECULAR_TRANSPORT

        self.is_ArakawaKonor = 0
        self.is_nonhydrostatic = 1
        self.is_compressible = 1

        self.compressibility = 0.0
        self.acoustic_timestep = 0
        self.acoustic_order = 0
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.gravity_strength = np.zeros((3))
        self.coriolis_strength = np.zeros((3))

        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        self.coriolis_strength[0] = self.omega * self.t_ref
        self.coriolis_strength[2] = self.omega * self.T_ref

        for i in range(3):
            if (self.gravity_strength[i] > np.finfo(np.float).eps) or (i == 1):
                self.i_gravity[i] = 1
                self.gravity_direction = i

            if (self.coriolis_strength[i] > np.finfo(np.float).eps):
                self.i_coriolis[i] = 1

        self.xmin = - 1.0 
        self.xmax =   1.0
        self.ymin =   0.0
        self.ymax =   1.0
        self.zmin = - 1.0
        self.zmax =   1.0

        self.u_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.WALL
        self.bdry_type[2] = BdryType.WALL

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL  = 0.5
        # self.CFL = 0.9
        # self.dtfixed0 = 1.9 / self.t_ref
        # self.dtfixed = 1.9 / self.t_ref

        # self.CFL = 1.0
        # self.dtfixed0 = 20.0 / self.t_ref
        # self.dtfixed = 20.0 / self.t_ref
        # self.dtfixed0 = 0.05
        # self.dtfixed = 0.05
        self.dtfixed0 = 100.0
        self.dtfixed = 100.0

        self.inx = 160+1
        self.iny = 80+1
        self.inz = 1

        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.eps_Machine = np.sqrt(np.finfo(np.float).eps)

        self.tout = np.arange(0.0,10.10,0.10)[10:]
        # self.tout = [10.0]

        self.stepmax = 10000

        self.blending_weight = 0./16
        self.blending_mean = 'rhoY' # 1.0, rhoY
        self.blending_conv = 'rho' #theta, rho

        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.initial_blending = False

        self.tol = 1.e-8
        self.max_iterations = 6000

        self.output_base_name = "_rising_bubble"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.tout[-1])
        
        aux = 'debug'
        self.aux = aux
        self.output_suffix = "_%i_%i_%.1f_%s" %(self.inx-1,self.iny-1,self.tout[-1],aux)

        # self.output_suffix += '_w=%i-%i' %(self.blending_weight*16.0,16.0-(self.blending_weight*16.0))

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_method

    def stratification_function(self, y):
        return 1.0

    @staticmethod
    def rhoe_method(rho,u,v,w,p,ud,th):
        Msq = ud.compressibility * ud.Msq

        gm1inv = th.gm1inv
        return (p * gm1inv + 0.5 * Msq * rho * (u*u + v*v + w*w))

def sol_init(Sol, mpv, elem, node, th, ud, seed=None):
    u0 = ud.u_wind_speed
    v0 = 0.0
    w0 = ud.w_wind_speed
    delth = 2.0         # [K]
    
    y0 = 0.2
    r0 = 0.2

    g = ud.gravity_strength[1]
    # print(ud.rho_ref)

    hydrostatic_state(mpv, elem, node, th, ud)

    x = elem.x
    y = elem.y

    x, y = np.meshgrid(x,y)

    if seed != None:
        np.random.seed(seed)
        # y0 += (np.random.random()-.5)/2.0
        # delth += 10.0*(np.random.random()-.5)
        delth += 10.0*(np.random.random())
    
    if 'truth' in ud.aux:
        np.random.seed(1234)
        # delth += 10.0*(np.random.random()-.5)
        delth += 10.0*(np.random.random())
    print(delth)
    
    r = np.sqrt((x)**2 + (y-y0)**2) / r0

    p = np.repeat(mpv.HydroState.p0.reshape(1,-1),elem.icx,axis=0)
    rhoY = np.repeat(mpv.HydroState.rhoY0.reshape(1,-1),elem.icx,axis=0)

    perturbation = (delth/300.0) * (np.cos(0.5 * np.pi * r)**2)
    perturbation[np.where(r > 1.0)] = 0.0
    rho = rhoY / (ud.stratification(y) + perturbation.T)

    x_idx = slice(None)
    y_idx = slice(None)

    u, v, w = u0, v0, w0

    Sol.rho[x_idx,y_idx] = rho
    Sol.rhou[x_idx,y_idx] = rho * u
    Sol.rhov[x_idx,y_idx] = rho * v
    Sol.rhow[x_idx,y_idx] = rho * w
    Sol.rhoe[x_idx,y_idx] = ud.rhoe(rho, u ,v, w, p, ud, th)
    Sol.rhoY[x_idx,y_idx] = rhoY

    p = mpv.HydroState_n.p0[0]
    rhoY = mpv.HydroState_n.rhoY0[0]
    mpv.p2_nodes[...] = (p - mpv.HydroState_n.p0[0]) / rhoY / ud.Msq

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    return Sol