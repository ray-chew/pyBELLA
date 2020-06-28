import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
from physics.hydrostatics import hydrostatic_state
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2
from physics.low_mach.second_projection import euler_backward_non_advective_impl_part

from scipy import signal

class UserData(object):
    NSPEC = 1

    grav = 0.0
    omega = 0.0

    R_gas = 287.4
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
        # self.Msq = 1.0

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

        self.xmin = - 5E+5
        self.xmax =   5E+5
        self.ymin = - 5E+5
        self.ymax =   5E+5
        self.zmin = - 0.5
        self.zmax =   0.5

        self.wind_speed = 0.0

        self.bdry_type_min = np.empty((3), dtype=object)
        self.bdry_type_max = np.empty((3), dtype=object)

        self.bdry_type_min[0] = BdryType.WALL
        self.bdry_type_min[1] = BdryType.WALL
        self.bdry_type_min[2] = BdryType.WALL
        self.bdry_type_max[0] = BdryType.WALL
        self.bdry_type_max[1] = BdryType.WALL
        self.bdry_type_max[2] = BdryType.WALL

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.WALL
        self.bdry_type[1] = BdryType.WALL
        self.bdry_type[2] = BdryType.WALL

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.95
        # self.CFL = 0.9 / 2.0
        self.dtfixed0 = 1000.0
        self.dtfixed = 1000.0

        self.inx = 150+1
        self.iny = 150+1
        self.inz = 1

        self.recovery_order = RecoveryOrder.SECOND
        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.kp = 0.0
        self.kz = 0.0
        self.km = 0.0
        self.kY = 0.0
        self.kZ = 0.0

        self.tol = 1.e-8
        self.max_iterations = 6000

        self.perturb_type = 'pos_perturb'
        self.blending_mean = 'rhoY' # 1.0, rhoY
        self.blending_conv = 'rho' #theta, rho

        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.blending_weight = 16./16

        self.initial_projection = True
        self.initial_impl_Euler = False

        self.column_preconditionr = False
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 0.0

        stepsize = 800
        # self.tout = np.arange(0,2E5+stepsize,stepsize)
        self.stepmax = 31
        self.tout = np.arange(0,self.stepmax+stepsize,stepsize)
        

        self.output_base_name = "_swe"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.tout[-1])
        
        # aux = 'posp_rloc'
        # aux += '_' + self.blending_conv + '_conv'
        # aux += '_' + self.blending_mean + '_mean'
        # aux = 'cb1_w=-6_debug'
        # self.output_suffix += '_w=%i-%i' %(self.blending_weight*16.0,16.0-(self.blending_weight*16.0))
        # aux = 'psinc_bal_debug'
        # self.output_suffix = "_%i_%i_%.1f_%s" %(self.inx-1,self.iny-1,self.tout[-1],aux)

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
    # wind speed
    u0 = 0.0
    v0 = 0.0
    w0 = 0.0

    # initial velocities
    u, v, w = 0.0, 0.0 , 0.0

    igs = elem.igs
    igy = igs[1]

    # g = ud.gravity_strength[1]
    g = 9.81    # [m/s^2]
    H = 100.0 * ud.h_ref  # [m]

    # dx = np.diff(elem.x).min()
    # dy = np.diff(elem.y).min()
    dx = 1E6/149
    dy = dx
    ud.tout = np.arange(0,500000.0,0.1*min(dx, dy)/np.sqrt(g*H)*20)
    # ud.dtfixed = 0.1*min(dx, dy)/np.sqrt(g*H) * 100.0

    i2 = (slice(igs[0],-igs[0]),slice(igs[1],-igs[1]))

    hydrostatic_state(mpv, elem, node, th, ud)

    # X,Y = np.meshgrid(elem.x[igs[0:-igs[0]]],elem.y[igs[1]:-igs[1]])
    X,Y = np.meshgrid(elem.x,elem.y)
    rho = np.ones_like(Sol.rho[i2]) * H
    # rho0 = np.copy(rho)
    # rho += np.exp(-((X-X.max()/2.7)**2/(2*(0.05E+6)**2) + (Y-Y.max()/4)**2/(2*(0.05E+6)**2)))
    rho += np.load('./output_swe/eta_list.npy')[0]
    # rhoH = rho + H

    # rhoH = rho + rho0

    Sol.rho[i2] = rho
    Sol.rhou[i2] = rho * u
    Sol.rhov[i2] = rho * v
    Sol.rhow[i2] = rho * w

    # Sol.rho[...] = rho
    # Sol.rhou[...] = rho * u
    # Sol.rhov[...] = rho * v
    # Sol.rhow[...] = rho * w

    if (ud.is_compressible) :
        Sol.rhoY[...] = 1.0
    else:
        Sol.rhoY[...] = 1.0

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    kernel = np.ones((2,2))
    kernel /= kernel.sum()
    rho_n = signal.convolve(Sol.rho, kernel, mode='valid')

    # print(rho_n.shape)
    # mpv.p2_nodes[igxn:-igxn,igyn:-igyn] = g *s rho_n**2 / 2.0
    mpv.p2_nodes[1:-1,1:-1] = g * rho_n**th.gamm / 2.0
    set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)
