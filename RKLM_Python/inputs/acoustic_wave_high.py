import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
# from discretization.time_discretization import SetTimeIntegratorParameters
# from physics.gas_dynamics.explicit import TimeIntegratorParams
from physics.gas_dynamics.eos import rhoe
from physics.hydrostatics import hydrostatic_state
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2, get_tau_y

class UserData(object):
    NSPEC = 1
    BOUY = 0

    grav = 0.0
    omega = 0.0*0.0001

    R_gas = 287.0
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
    T_ref = 353.048780488
    p_ref = 101325
    u_ref = h_ref / t_ref
    rho_ref = p_ref / (R_gas * T_ref)

    Nsq_ref = 0.0e-4

    scale_factor = 1.0

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
        self.viscm = self.viscm * self.t_ref / (self.h_ref * self.h_ref)
        self.viscbm = self.viscbm * self.t_ref / (self.h_ref * self.h_ref)
        self.visct = self.visct * self.t_ref / (self.h_ref * self.h_ref)
        self.viscbt = self.viscbt * self.t_ref / (self.h_ref * self.h_ref)
        self.cond = self.cond * self.t_ref / (self.h_ref * self.h_ref * self.R_gas)

        self.is_nonhydrostatic = 1
        self.is_compressible = 1
        self.is_ArakawaKonor = 0
        self.compressibility = 0.0
        self.acoustic_timestep = 0
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

        self.xmin =   0.0 * self.scale_factor
        self.xmax =   1.0 * self.scale_factor
        self.ymin =   0.0 * self.scale_factor
        self.ymax =   1.0 * self.scale_factor
        self.zmin = - 1.0
        self.zmax =   1.0

        self.u_wind_speed = 0.0 / self.u_ref
        self.v_wind_speed = 0.0 / self.u_ref
        self.w_wind_speed = 0.0 / self.u_ref

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

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.RAYLEIGH
        self.bdry_type[2] = BdryType.WALL

        ##########################################
        # NUMERICS
        ##########################################
        # self.CFL  = 77.0
        self.CFL = 0.45
        self.dtfixed0 = 0.0000668205
        self.dtfixed = 0.0000668205

        # self.inx = 64+1
        # self.iny = 32+1
        self.inx = 64+1
        self.iny = 64+1
        self.inz = 1

        if self.bdry_type[1].value == 'radiation':
            self.inbcy = self.iny - 1
            self.iny += self.inbcy

            # tentative workaround
            self.bcy = self.ymax
            self.ymax *= 2.0        

        self.recovery_order = RecoveryOrder.SECOND
        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.tol = 1.e-8
        self.max_iterations = 1000

        self.blending_weight = 0./16
        self.blending_mean = 'rhoY' # 1.0, rhoY
        self.blending_conv = 'rho' # theta, rho
        self.blending_type = 'half'
        
        self.continuous_blending = False
        self.initial_blending = False
        self.no_of_pi_initial = 0
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.initial_projection = False

        # self.tout[0] =  0.00267282
        self.tout = [0.00267282]
        # self.tout = np.arange(0.0,0.00267282,2.5e-5)
        # self.tout = np.arange(0.0,0.002+2.5e-5,2.5e-5)[1:]
        # self.tout[1] = -1.0

        self.stepmax = 999

        self.output_timesteps = True

        self.output_base_name = "_acoustic_wave"
        self.output_name_psinc = "_low_mach_gravity_psinc"
        self.output_name_comp = "_low_mach_gravity_comp"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%.5f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%.5f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%.5f" %(self.inx-1,self.iny-1,self.tout[-1])

        self.aux = 'test'
        
        self.stratification = self.stratification_function

    def stratification_function(self, y):
        if type(y) == float:
            return 1.0
        else:
            return np.ones((y.shape))

def sol_init(Sol, mpv, elem, node, th, ud):
    u0 = ud.u_wind_speed
    v0 = 0.0
    w0 = 0.0
    del0 = 0.05
    wn = 2.0 * np.pi / (ud.xmax - ud.xmin)

    Ma = np.sqrt(ud.Msq)

    hydrostatic_state(mpv, elem, node, th, ud)

    if ud.bdry_type[1].value == 'radiation':
        tcy, tcn = get_tau_y(ud, elem, node, 2.5)

    x_idx = slice(None)
    x = elem.x[x_idx].reshape(-1,1)
    y_idx = slice(elem.igy,-elem.igy+1)

    v, w = v0, w0
    p = mpv.HydroState.p0[y_idx] * (1.0 + del0 * np.sin(wn * x))**(2.0 * th.gamm * th.gm1inv)
    rhoY = p**(th.gamminv)
    rho = np.copy(rhoY)

    c = np.sqrt(th.gamm * np.divide(p , rho))
    u = u0 + (p - mpv.HydroState.p0[y_idx]) / (rho * c) / Ma
    
    # u_shape = u.shape
    # u = np.cumsum(u.T.reshape(-1)).reshape(u_shape)

    Sol.rho[x_idx,y_idx] = rho
    Sol.rhou[x_idx,y_idx] = rho * u
    Sol.rhov[x_idx,y_idx] = rho * v
    Sol.rhow[x_idx,y_idx] = rho * w
    Sol.rhoe[x_idx,y_idx] = rhoe(rho, u ,v, w, p, ud, th)
    Sol.rhoY[x_idx,y_idx] = rhoY

    mpv.p2_cells[x_idx,y_idx] = (p**th.Gamma - 1.0) / ud.Msq

    Sol.rhoX[x_idx,y_idx] = Sol.rho[x_idx,y_idx] * (Sol.rho[x_idx,y_idx] / Sol.rhoY[x_idx,y_idx] - mpv.HydroState.S0[y_idx])

    x = node.x[x_idx][:-1].reshape(-1,1)
    p = mpv.HydroState_n.p0[y_idx] * (1.0 + del0 * np.sin(wn * x))**(2.0 * th.gamm * th.gm1inv)
    mpv.p2_nodes[:-1,y_idx] = (p**th.Gamma - 1.0) / ud.Msq

    ud.initial_projection = False

    set_ghostcells_p2(mpv.p2_cells,elem,ud)
    set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    ud.is_nonhydrostasy = 1.0
    ud.compressibility = 1.0 if ud.is_compressible == 1 else 0.0

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)
