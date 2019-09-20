import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
from discretization.time_discretization import SetTimeIntegratorParameters
from physics.gas_dynamics.explicit import TimeIntegratorParams
from physics.hydrostatics import hydrostatic_state, hydrostatic_column, hydrostatic_initial_pressure
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2
from management.variable import States

class UserData(object):
    NSPEC = 1
    BOUY = 0

    grav = 10.0
    omega = 0.0*0.0001

    R_gas = 287.4
    R_vap = 461.0
    Q_vap = 2.53e+06
    gamma = 1.4

    viscm = 0.0
    viscbm = 0.0
    visct = 0.0
    viscbt = 0.0
    cond = 0.0

    h_ref = 10000.0
    t_ref = 100.0
    T_ref = 300.00
    p_ref = 8.61 * 1e4
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
        self.viscm = self.viscm * self.t_ref / (self.h_ref * self.h_ref)
        self.viscbm = self.viscbm * self.t_ref / (self.h_ref * self.h_ref)
        self.visct = self.visct * self.t_ref / (self.h_ref * self.h_ref)
        self.viscbt = self.viscbt * self.t_ref / (self.h_ref * self.h_ref)
        self.cond = self.cond * self.t_ref / (self.h_ref * self.h_ref * self.R_gas)

        self.is_nonhydrostatic = 1
        self.is_compressible = 1
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

        self.xmin = - 1.0 
        self.xmax =   1.0
        self.ymin =   0.0
        self.ymax =   1.0
        self.zmin = - 1.0
        self.zmax =   1.0

        self.wind_speed = 0.0
        self.wind_shear = -0.0
        self.hill_shape = HillShapes.AGNESI
        self.hill_height = 0.0 * 0.096447
        self.hill_length_scale = 0.1535

        self.bdry_type_min = np.empty((3), dtype=object)
        self.bdry_type_max = np.empty((3), dtype=object)

        self.bdry_type_min[0] = BdryType.PERIODIC
        self.bdry_type_min[1] = BdryType.WALL
        self.bdry_type_min[2] = BdryType.WALL
        self.bdry_type_max[0] = BdryType.PERIODIC
        self.bdry_type_max[1] = BdryType.WALL
        self.bdry_type_max[2] = BdryType.WALL

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.WALL
        self.bdry_type[2] = BdryType.WALL

        self.absorber = False # 0 == WRONG == FALSE 
        self.bottom_theta_bc = BottomBC.BOTTOM_BC_DEFAULT

        ##########################################
        # NUMERICS
        ##########################################

        self.time_integrator = TimeIntegrator.SI_MIDPT
        self.advec_time_integrator = TimeIntegrator.STRANG
        self.CFL  = 1.0
        self.dtfixed0 = 17.0 / self.t_ref
        self.dtfixed = 17.0 / self.t_ref

        self.tips = TimeIntegratorParams()
        SetTimeIntegratorParameters(self)

        self.inx = 160+1
        self.iny = 80+1
        self.inz = 1

        self.recovery_order = RecoveryOrder.SECOND
        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.kp = 1.4
        self.kz = 1.4
        self.km = 1.4
        self.kY = 1.4
        self.kZ = 1.4

        self.ncache = 154

        tol = 1.e-11 if self.is_compressible == 0 else 1.e-11 * 0.01

        self.flux_correction_precision = tol
        self.flux_correction_local_precision = tol
        self.second_projection_precision = tol
        self.second_projection_local_precision = tol
        self.flux_correction_max_iterations = 6000
        self.second_projection_max_iterations = 6000

        self.initial_projection = False
        self.initial_impl_Euler = False

        self.column_preconditionr = True
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 1.0

        self.eps_Machine = np.sqrt(np.finfo(np.float).eps)

        self.tout =  1000.0 / self.t_ref
        # self.tout[0] =  self.scale_factor * 1.0 * 3000.0 / self.t_ref
        # self.tout[1] = -1.0

        # self.stepmax = 10000
        self.stepmax = 1000

        self.write_stdout = True
        self.write_stdout_period = 1
        self.write_file = True
        self.write_file_period = 30
        self.file_format = 'HDF'

        self.output_base_name = "_rising_bubble"
        self.output_name_psinc = "_low_mach_gravity_psinc"
        self.output_name_comp = "_low_mach_gravity_comp"

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_method

    def stratification_function(self, y):
        return 1.0

    @staticmethod
    def rhoe_method(rho,u,v,w,p,ud,th):
        Msq = ud.compressibility * ud.Msq

        gm1inv = th.gm1inv
        return (p * gm1inv + 0.5 * Msq * rho * (u*u + v*v + w*w))

def sol_init(Sol, mpv, elem, node, th, ud):
    u0 = ud.wind_speed
    v0 = 0.0
    w0 = 0.0
    delth = 2.0
    
    y0 = 0.2
    r0 = 0.2

    g = ud.gravity_strength[1]

    hydrostatic_state(mpv, elem, node, th, ud)

    x = elem.x.reshape(-1,1)
    y = elem.y.reshape(1,-1)

    r = np.sqrt(x**2 + (y-y0)**2 ) / r0
    r = r.reshape(Sol.rho.shape)
    p = np.repeat(mpv.HydroState.p0[0].reshape(1,-1),elem.icx,axis=0)
    rhoY = np.repeat(mpv.HydroState.rhoY0[0].reshape(1,-1),elem.icx,axis=0)
    
    perturbation = ((delth/300.) * np.cos(0.5 * np.pi * r)**2)
    perturbation[np.where(r >= 1.0)] = 0.0
    rho = rhoY / ud.stratification(y) - perturbation

    # Y = ud.stratification(y) + delth * ud.molly(x) * np.sin(np.pi * y) / (1.0 + (x-xc)**2 / (a**2))

    # xn = node.x[:-1].reshape(-1,1)
    # yn = node.y[:-1].reshape(1,-1)

    # Yn = ud.stratification(yn) + delth * ud.molly(xn) * np.sin(np.pi * yn) / (1.0 + (xn-xc)**2 / (a**2))

    # hydrostatic_column(HySt, HyStn, Y, Yn, elem, node, th, ud)
    # print('mpv.HydroState_n.Y0 = ', mpv.HydroState_n.Y0[0])
    # print('HyStn = ', HyStn.Y0[0])
    x_idx = slice(None)
    y_idx = slice(None)
    # xc_idx = slice(0,-1)
    # yc_idx = slice(0,-1)
    # c_idx = (xc_idx, yc_idx)

    u, v, w = u0, v0, w0
    # if ud.is_compressible:
    #     p = HySt.p0[:,y_idx][c_idx]
    #     rhoY = HySt.rhoY0[:,y_idx][c_idx]
    # else:
    #     p = mpv.HydroState.p0[:, y_idx]
    #     rhoY = mpv.HydroState.rhoY0[:, y_idx]

    # rho = rhoY / Y[:,y_idx]
    Sol.rho[x_idx,y_idx] = rho
    Sol.rhou[x_idx,y_idx] = rho * u
    Sol.rhov[x_idx,y_idx] = rho * v
    Sol.rhow[x_idx,y_idx] = rho * w
    Sol.rhoe[x_idx,y_idx] = ud.rhoe(rho, u ,v, w, p, ud, th)
    Sol.rhoY[x_idx,y_idx] = rhoY

    # p = mpv.HydroState_n.p0[0]
    # rhoY = mpv.HydroState_n.rhoY0[0]
    # mpv.p2_nodes = (p - mpv.HydroState_n.p0[0]) / rhoY / ud.Msq

    # mpv.p2_cells[x_idx,y_idx] = HySt.p20[x_idx,y_idx][c_idx]
    # mpv.p2_cells[x_idx,y_idx] -= mpv.HydroState.p20[0,y_idx]
    # print(mpv.p2_cells[x_idx,y_idx] - mpv.HydroState.p20[0,y_idx])
    # print(mpv.p2_cells[0])
    # Sol.rhoX[x_idx,y_idx] = Sol.rho[x_idx,y_idx] * (1.0 / Y[:, y_idx] - mpv.HydroState.S0[0, y_idx])

    # mpv.p2_nodes[:,elem.igy:-elem.igy] = HyStn.p20[:,elem.igy:-elem.igy]

    # hydrostatic_initial_pressure(Sol,mpv,elem,node,ud,th)
    # print(mpv.p2_nodes[103,:10])
    ud.is_nonhydrostasy = 0.0
    ud.compressibility = 1.0

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    return Sol