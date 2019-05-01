import numpy as np
from input.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
from numerics_fundamentals.discretization.time_discretization import SetTimeIntegratorParameters
from physics.gas_dynamics.explicit import TimeIntegratorParams

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

    h_ref = 100
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

        self.mol_trans = MolecularTransport.NO_MOLECULAR_TRANSPORT
        self.viscm = self.viscm * self.t_ref / (self.h_ref * self.h_ref)
        self.viscbm = self.viscbm * self.t_ref / (self.h_ref * self.h_ref)
        self.visct = self.visct * self.t_ref / (self.h_ref * self.h_ref)
        self.viscbt = self.viscbt * self.t_ref / (self.h_ref * self.h_ref)
        self.cond = self.cond * self.t_ref / (self.h_ref * self.h_ref * self.R_gas)

        self.is_nonhydrostatic = 1
        self.is_compressible = 1
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
        self.bdry_type_min[2] = BdryType.WALL
        self.bdry_type_max[0] = BdryType.PERIODIC
        self.bdry_type_max[1] = BdryType.PERIODIC
        self.bdry_type_max[2] = BdryType.WALL

        self.absorber = 0 # 0 == WRONG == FALSE 
        self.bottom_theta_bc = BottomBC.BOTTOM_BC_DEFAULT

        ##########################################
        # NUMERICS
        ##########################################

        self.time_integrator = TimeIntegrator.SI_MIDPT
        self.advec_time_integrator = TimeIntegrator.STRANG
        self.CFL  = 0.9/2.0
        self.dtfixed0 = 2.1 * 1.200930e-2
        self.dtfixed = 2.1 * 1.200930e-2

        self.tips = TimeIntegratorParams()
        SetTimeIntegratorParameters(self)

        self.inx = 48+1
        self.iny = 48+1
        self.inz = 1

        self.recovery_order = RecoveryOrder.SECOND
        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.kp = 0.0
        self.kz = 0.0
        self.km = 0.0
        self.kY = 0.0
        self.kZ = 0.0

        self.ncache = 201

        tol = 1.e-10

        self.flux_correction_precision = tol
        self.flux_correction_local_precision = tol
        self.second_projection_precision = tol
        self.second_projection_local_precision = tol
        self.flux_correction_max_iterations = 6000
        self.second_projection_max_iterations = 6000

        self.initial_projection = True
        self.initial_impl_Euler = False

        self.column_preconditionr = False
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 0.0

        self.eps_Machine = np.sqrt(np.finfo(np.float).eps)

        self.tout[0] =  1.0
        self.tout[1] = -1.0

        self.stepmax = 20000

        self.write_stdout = True
        self.write_stdout_period = 1
        self.write_file = True
        self.write_file_period = 10000
        self.file_format = 'HDF'

        self.n_time_series = 500

        self.output_base_folder = "output/"
        self.output_folder_name_psinc = "low_Mach_gravity_psinc/"
        self.output_folder_name_comp = "low_Mach_gravity_comp"



        



