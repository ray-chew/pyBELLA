import numpy as np
import inputs.enum_bdry as enum_bdry
import management.enumerator as enumerator
import physics.hydrostatics as hydrostatics
import inputs.boundary as boundary
import management.variable as variable

class UserData(object):
    NSPEC = 1

    grav = 9.81
    omega = 1.0 * 0.0001
    R_gas = 287.4
    R_vap = 461.00
    Q_vap = 2.53e06
    gamma = 1.4

    h_ref = 1000
    t_ref = 100
    T_ref = 300.00
    p_ref = 1e5
    u_ref = h_ref / t_ref
    rho_ref = p_ref / (R_gas * T_ref)

    Nsq_ref = 0.00011772

    i_gravity = np.zeros((3))
    i_coriolis = np.zeros((3))

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

        self.is_ArakawaKonor = 0
        self.is_nonhydrostatic = 1
        self.is_compressible = 1
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
                self.gravity_direction = 1

            if (self.coriolis_strength[i] > np.finfo(np.float).eps):
                self.i_coriolis[i] = 1

        self.xmin =   0.0e6 / self.h_ref
        self.xmax =  10.0e6 / self.h_ref
        self.ymin =   0.0e0 / self.h_ref
        self.ymax =  18.0e3 / self.h_ref
        self.zmin = - 8.0e6 / self.h_ref
        self.zmax =   8.0e6 / self.h_ref

        self.wind_speed = 0.0 * 20.0 / self.u_ref
        self.wind_shear = -0.0
        self.hill_shape = enumerator.HillShapes.AGNESI
        self.hill_height = 0.0 * 0.096447
        self.hill_length_scale = 0.1535

        self.bdry_type_min = np.empty((3), dtype=object)
        self.bdry_type_max = np.empty((3), dtype=object)

        self.bdry_type_min[0] = enum_bdry.BdryType.PERIODIC
        self.bdry_type_min[1] = enum_bdry.BdryType.WALL
        self.bdry_type_min[2] = enum_bdry.BdryType.PERIODIC
        self.bdry_type_max[0] = enum_bdry.BdryType.PERIODIC
        self.bdry_type_max[1] = enum_bdry.BdryType.WALL
        self.bdry_type_max[2] = enum_bdry.BdryType.PERIODIC

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = enum_bdry.BdryType.PERIODIC
        self.bdry_type[1] = enum_bdry.BdryType.WALL
        self.bdry_type[2] = enum_bdry.BdryType.PERIODIC

        self.absorber = 0 # 0 == WRONG == FALSE 
        self.bottom_theta_bc = enumerator.BottomBC.BOTTOM_BC_DEFAULT

        ##########################################
        # NUMERICS
        ##########################################
        self.time_integrator = enumerator.TimeIntegrator.SI_MIDPT
        self.advec_time_integrator = enumerator.TimeIntegrator.STRANG
        self.CFL  = 0.4
        self.dtfixed0 = 8.0
        self.dtfixed = 8.0

        # self.tips = TimeIntegratorParams()
        # SetTimeIntegratorParameters(self)

        self.inx = 65+1
        self.iny = 33+1
        self.inz = 65+1

        self.recovery_order = enumerator.RecoveryOrder.SECOND
        self.limiter_type_scalars = enumerator.LimiterType.NONE
        self.limiter_type_velocity = enumerator.LimiterType.NONE

        self.kp = 1.4
        self.kz = 1.4
        self.km = 1.4
        self.kY = 1.4
        self.kZ = 1.4

        self.continuous_blending = False
        self.no_of_pi_initial = 0
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        tol = 1.e-8
        self.flux_correction_precision = tol
        self.flux_correction_local_precision = tol
        self.second_projection_precision = tol
        self.second_projection_local_precision = tol
        self.flux_correction_max_iterations = 6000
        self.second_projection_max_iterations = 6000
        self.initial_projection = False

        self.column_preconditionr = 1
        # self.synchronize_nodal_pressure = False
        # self.synchronize_weight = 0.0

        self.eps_Machine = np.sqrt(np.finfo(np.float).eps)

        self.tout =  [2.5e1]
        # self.tout[1] = -1.0

        self.stepmax = 20000

        self.output_base_name = "_baroclinic_instability_periodic"
        self.output_name_psinc = "_low_mach_gravity_psinc"
        self.output_name_comp = "_low_mach_gravity_comp"

    def sol_init(Sol, mpv, elem, node, th, ud, seed=None):
        None