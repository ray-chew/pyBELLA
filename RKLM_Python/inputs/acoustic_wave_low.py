# import numpy as np
# from inputs.enum_bdry import BdryType
# from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
# from numerics_fundamentals.discretization.time_discretization import SetTimeIntegratorParameters
# from physics.gas_dynamics.explicit import TimeIntegratorParams
# from physics.hydrostatics.hydrostatics import HydrostaticStates
# from inputs.boundary import set_wall_rhoYflux, set_explicit_boundary_data, set_ghostcells_p2

# class UserData(object):
#     NSPEC = 1
#     BOUY = 0

#     grav = 0.0
#     omega = 0.0*0.0001

#     R_gas = 287.4
#     R_vap = 461.0
#     Q_vap = 2.53e+06
#     gamma = 2.0

#     viscm = 0.0
#     viscbm = 0.0
#     visct = 0.0
#     viscbt = 0.0
#     cond = 0.0

#     h_ref = 1.0
#     t_ref = 1.0
#     T_ref = 353.048780488
#     p_ref = 101325
#     u_ref = h_ref / t_ref
#     rho_ref = p_ref / (R_gas * T_ref)

#     Nsq_ref = 0.0e-4

#     scale_factor = 1.0

#     i_gravity = np.zeros((3))
#     i_coriolis = np.zeros((3))

#     tout = np.zeros((2))

#     def __init__(self):
#         self.h_ref = self.h_ref
#         self.t_ref = self.t_ref
#         self.T_ref = self.T_ref
#         self.p_ref = self.p_ref
#         self.rho_ref = self.rho_ref
#         self.u_ref = self.u_ref
#         self.Nsq_ref = self.Nsq_ref
#         self.g_ref = self.grav
#         self.gamm = self.gamma
#         self.Rg_over_Rv = self.R_gas / self.R_vap
#         self.Q = self.Q_vap / (self.R_gas * self.T_ref)

#         self.nspec = self.NSPEC

#         self.mol_trans = MolecularTransport.NO_MOLECULAR_TRANSPORT
#         self.viscm = self.viscm * self.t_ref / (self.h_ref * self.h_ref)
#         self.viscbm = self.viscbm * self.t_ref / (self.h_ref * self.h_ref)
#         self.visct = self.visct * self.t_ref / (self.h_ref * self.h_ref)
#         self.viscbt = self.viscbt * self.t_ref / (self.h_ref * self.h_ref)
#         self.cond = self.cond * self.t_ref / (self.h_ref * self.h_ref * self.R_gas)

#         self.is_nonhydrostatic = 1
#         self.is_compressible = 1
#         self.acoustic_timestep = 1
#         self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

#         self.gravity_strength = np.zeros((3))
#         self.coriolis_strength = np.zeros((3))

#         self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
#         self.coriolis_strength[0] = self.omega * self.t_ref
#         self.coriolis_strength[2] = self.omega * self.T_ref

#         for i in range(3):
#             if (self.gravity_strength[i] > np.finfo(np.float).eps) or (i == 1):
#                 self.i_gravity[i] = 1
#                 self.gravity_direction = 1

#             if (self.coriolis_strength[i] > np.finfo(np.float).eps):
#                 self.i_coriolis[i] = 1

#         self.xmin =   0.0 * self.scale_factor
#         self.xmax =   1.0 * self.scale_factor
#         self.ymin =   0.0 * self.scale_factor
#         self.ymax =   1.0 * self.scale_factor
#         self.zmin = - 1.0 * self.scale_factor
#         self.zmax =   1.0 * self.scale_factor

#         self.wind_speed = 1.0 / self.u_ref
#         self.wind_shear = -0.0
#         self.hill_shape = HillShapes.AGNESI
#         self.hill_height = 0.0 * 0.096447
#         self.hill_length_scale = 0.1535

#         # self.xmin = - 5000. / self.h_ref
#         # self.xmax =   5000. / self.h_ref
#         # self.ymin = - 5000. / self.h_ref / 8.0
#         # self.ymax =   5000. / self.h_ref / 8.0
#         # self.zmin = - 5000. / self.h_ref
#         # self.zmax =   5000. / self.h_ref

#         # self.wind_speed = 1.0 * 10.0 / self.u_ref
#         # self.wind_shear = -0.0
#         # self.hill_height = 0.0
#         # self.hill_length_scale = 99999.9

#         self.bdry_type_min = np.empty((3), dtype=object)
#         self.bdry_type_max = np.empty((3), dtype=object)

#         self.bdry_type_min[0] = BdryType.PERIODIC
#         self.bdry_type_min[1] = BdryType.WALL
#         self.bdry_type_min[2] = BdryType.WALL
#         self.bdry_type_max[0] = BdryType.PERIODIC
#         self.bdry_type_max[1] = BdryType.WALL
#         self.bdry_type_max[2] = BdryType.WALL

#         self.bdry_type = np.empty((3), dtype=object)
#         self.bdry_type[0] = BdryType.PERIODIC
#         self.bdry_type[1] = BdryType.WALL
#         self.bdry_type[2] = BdryType.WALL

#         self.absorber = 0 # 0 == WRONG == FALSE 
#         self.bottom_theta_bc = BottomBC.BOTTOM_BC_DEFAULT

#         ##########################################
#         # NUMERICS
#         ##########################################

#         self.time_integrator = TimeIntegrator.SI_MIDPT
#         self.advec_time_integrator = TimeIntegrator.STRANG
#         self.CFL  = 10.0
#         self.dtfixed0 = 0.0000668205
#         self.dtfixed = 0.0000668205

#         self.tips = TimeIntegratorParams()
#         SetTimeIntegratorParameters(self)

#         self.inx = 256+1
#         self.iny = 10+1
#         self.inz = 1

#         self.recovery_order = RecoveryOrder.SECOND
#         self.limiter_type_scalars = LimiterType.NONE
#         self.limiter_type_velocity = LimiterType.NONE

#         self.kp = 1.4
#         self.kz = 1.4
#         self.km = 1.4
#         self.kY = 1.4
#         self.kZ = 1.4

#         self.ncache = 175

#         tol = 1.e-8

#         self.flux_correction_precision = tol
#         self.flux_correction_local_precision = tol
#         self.second_projection_precision = tol
#         self.second_projection_local_precision = tol
#         self.flux_correction_max_iterations = 6000
#         self.second_projection_max_iterations = 6000

#         self.initial_projection = False
#         self.initial_impl_Euler = False

#         self.column_preconditionr = True
#         self.synchronize_nodal_pressure = False
#         self.synchronize_weight = 0.0

#         self.eps_Machine = np.sqrt(np.finfo(np.float).eps)

#         self.tout[0] =  0.00267282
#         self.tout[1] = -1.0

#         self.stepmax = 40

#         self.write_stdout = True
#         self.write_stdout_period = 1
#         self.write_file = True
#         self.write_file_period = 40
#         self.file_format = 'HDF'

#         self.n_time_series = 500

#         self.output_base_folder = "output/"
#         self.output_folder_name_psinc = "low_Mach_gravity_psinc/"
#         self.output_folder_name_comp = "low_Mach_gravity_comp"

#         self.stratification = self.stratification_function

#     def stratification_function(self, y):
#         return 1.0

# def sol_init(Sol, Sol0, mpv, bdry, elem, node, th, ud):
#     u0 = ud.wind_speed
#     v0 = 0.0
#     w0 = 0.0
#     del0 = 0.01
#     xc = 0.0
#     a = 0.125 * (ud.xmax - ud.xmin)
#     wn = 2.0 * np.pi / (ud.xmax - ud.xmin)

#     Ma = np.sqrt(ud.Msq)

#     HydrostaticStates(mpv, elem, node, th, ud)

#     x_idx = slice(None)
#     x = elem.x
#     y_idx = slice(elem.igy,-elem.igy)
#     y = elem.y

#     u,v,w = u0,v0,w0

#     p = mpv.HydroState.p0 * (1.0 + del0 * np.sin(wn * x))**(2.0 * th.gamm * th.gamm1inv)
#     rhoY = p**(th.gamminv)
#     rho = np.copy(rhoY)
#     c = np.sqrt(th.gamm * p / rho)

#     u = (p - mpv.HydroState.p0) / (rho * c) / Ma
#     u_shape = u.shape
#     u = np.cumsum(u.T.reshape(-1)).reshape(u_shape).T

#     Sol.rho[y_idx,x_idx] = rho
#     Sol.rhou[y_idx,x_idx] = rho * u
#     Sol.rhov[y_idx,x_idx] = rho * v
#     Sol.rhow[y_idx,x_idx] = rho * w
#     Sol.rhoe[y_idx,x_idx] = rhoe(rho, u ,v, w, p, ud, th)

#     mpv.p2_cells = (p**th.Gamma - 1.0) / ud.Msq


# def T_from_p_rho(p, rho):
#     return np.divide(p,rho)

# def rhoe(rho,u,v,w,p,ud,th):
#     Msq = ud.compressibility * ud.Msq
#     gm1inv = th.gm1inv

#     return p * gm1inv + 0.5 * Msq * rho * (u**2 + v**2 + w**2)