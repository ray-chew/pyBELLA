import numpy as np

from dycore.utils.options import BdryType, LimiterType
from sim_params import global_constants


class UserDataInit(object):
    """
    Loads user defined initial conditions. Specifically, all attributes of the class object defined in the initial condition is overwritten.

    Attributes
    ----------
    **kwargs: class object

    """

    def __init__(self,**kwargs):
        for key, value in global_constants.items():
            setattr(self, key, value)

        if len(kwargs) > 0:
            for key, value in kwargs.items():
                setattr(self, key, value)

        ##########################################
        # SPATIAL GRID
        ##########################################
        self.inx = 64+1
        self.iny = 64+1
        self.inz = 1

        self.xmin = - 1.0 
        self.xmax =   1.0
        self.ymin =   0.0
        self.ymax =   1.0
        self.zmin = - 1.0
        self.zmax =   1.0

        ##########################################
        # BOUNDARY CONDITIONS
        ##########################################
        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.WALL
        self.bdry_type[2] = BdryType.WALL

        ##########################################
        # TEMPORAL
        ##########################################
        self.CFL  = 0.5
        self.dtfixed0 = 100.0
        self.dtfixed = 100.0

        self.acoustic_timestep = 0

        self.tout = np.arange(0.0,1.01,0.01)[10:]
        self.stepmax = 10000

        ##########################################
        # MODEL REGIMES
        ##########################################
        self.is_ArakawaKonor = 0
        self.is_nonhydrostatic = 1
        self.is_compressible = 1

        self.compressibility = 0.0

        ##########################################
        # BACKGROUND WIND
        ##########################################
        self.u_wind_speed = 0.0
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        ##########################################
        # NUMERICS
        ##########################################   
        # Advection limiter types     
        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        # Iterative solver
        self.tol = 1.e-8
        self.max_iterations = 6000

        self.stratification = self.stratification_function

        ##########################################
        # BLENDING
        ##########################################
        self.blending_weight = 0./16
        self.blending_mean = 'rhoY' # 1.0, rhoY
        self.blending_conv = 'rho' # theta, rho
        self.blending_type = 'half'

        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.initial_blending = False

    def compute_u_ref(self):
        self.u_ref = self.h_ref / self.t_ref

    def compute_Msq(self):
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

    def compute_gravity(self):
        self.gravity_strength = np.zeros((3))
        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)

        for i in range(3):
            if (self.gravity_strength[i] > 0.0) or (i == 1):
                self.i_gravity[i] = 1
                self.gravity_direction = i

    def compute_coriolis(self):
        self.coriolis_strength = np.zeros((3))
        self.coriolis_strength[0] = self.omega * self.t_ref
        self.coriolis_strength[2] = self.omega * self.t_ref


    def compute_cp_gas(self):
        self.cp_gas = self.gamma * self.R_gas / (self.gamma-1.0)
        self.compute_N_ref()

    def compute_rho_ref(self):
        self.rho_ref = self.p_ref / (self.R_gas * self.T_ref)

    def compute_N_ref(self):
        self.N_ref = self.grav / np.sqrt(self.cp_gas * self.T_ref)
        self.Nsq_ref = self.N_ref * self.N_ref

    def compute_Cs(self):
        self.Cs = np.sqrt(self.gamma * self.R_gas * self.T_ref)
    

    @staticmethod
    def stratification_function():
        return 1.0


    def update_ud(self, obj):
        for key, value in obj.items():
            setattr(self, key, value)

    ##########################################
    # SETTER FUNCTIONS
    ##########################################

    # Msq arguments
    @property
    def u_ref(self):
        return self._u_ref
    
    @u_ref.setter
    def u_ref(self, val):
        self._u_ref = val
        self.compute_gravity()

    # gravity and Msq arguments
    @property
    def R_gas(self):
        return self._R_gas
    
    @R_gas.setter
    def R_gas(self, val):
        self._R_gas = val
        self.compute_gravity()
        self.compute_Msq()

        self.compute_cp_gas()
        self.compute_rho_ref()
        self.compute_Cs()

    @property
    def T_ref(self):
        return self._T_ref
    
    @T_ref.setter
    def T_ref(self, val):
        self._T_ref = val
        self.compute_gravity()
        self.compute_Msq()

        self.compute_rho_ref()
        self.compute_N_ref()
        self.compute_Cs()

    # gravity arguments
    def grav(self):
        return self._grav
    
    @grav.setter
    def grav(self, val):
        self._grav = val
        self.compute_gravity()
        self.compute_N_ref()

    @property
    def h_ref(self):
        return self._h_ref
    
    @h_ref.setter
    def h_ref(self, val):
        self._h_ref = val
        self.compute_u_ref()
        self.compute_gravity()

    # coriolis arguments
    @property
    def t_ref(self):
        return self._t_ref
    
    @t_ref.setter
    def t_ref(self, val):
        self._t_ref = val
        self.compute_u_ref()
        self.compute_coriolis()

    @property
    def omega(self):
        return self._omega
    
    @omega.setter
    def omega(self, val):
        self._omega = val
        self.compute_coriolis()


    # Cs and cp_gas argument
    @property
    def gamma(self):
        return self._gamma
    
    @gamma.setter
    def gamma(self, val):
        self._gamma = val
        self.compute_cp_gas()
        self.compute_Cs()

    # rho_ref argument
    @property
    def rho_ref(self):
        return self._rho_ref
    
    @rho_ref.setter
    def rho_ref(self, val):
        self._rho_ref = val
        self.compute_rho_ref()