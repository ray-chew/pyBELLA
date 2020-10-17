import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
from physics.hydrostatics import hydrostatic_state
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2
from physics.low_mach.second_projection import euler_backward_non_advective_impl_part

from scipy.interpolate import RectBivariateSpline
from scipy import signal

class UserData(object):
    NSPEC = 1

    grav = 0.0
    omega = 6.147 * 1E-5 # 360 # (24.0 * 60.0**2)    # [rad/s]
    # omega = 0.0001

    R_gas = 1.0
    R_vap = 461.0
    Q_vap = 2.53e+06
    gamma = 2.0

    viscm = 0.0
    viscbm = 0.0
    visct = 0.0
    viscbt = 0.0
    cond = 0.0

    h_ref = 1000.0        # [m]
    t_ref = 1200.0        # [day] -> [s]
    T_ref = 1.0
    p_ref = 1.0

    R_vap = 1.0
    Q_vap = 1.0

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
        self.R_gas = self.R_gas

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

        self.xmin =   0.0
        self.xmax =   5000.0*1E3 / self.h_ref
        self.ymin = - 0.5 * self.h_ref
        self.ymax =   0.5 * self.h_ref
        self.zmin =   0.0
        self.zmax =   4330.0*1E3 / self.h_ref

        self.u_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type_min = np.empty((3), dtype=object)
        self.bdry_type_max = np.empty((3), dtype=object)

        self.bdry_type_min[0] = BdryType.PERIODIC
        self.bdry_type_min[1] = BdryType.PERIODIC
        self.bdry_type_min[2] = BdryType.PERIODIC
        self.bdry_type_max[0] = BdryType.PERIODIC
        self.bdry_type_max[1] = BdryType.PERIODIC
        self.bdry_type_max[2] = BdryType.PERIODIC

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.PERIODIC
        self.bdry_type[2] = BdryType.PERIODIC

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL = 0.45
        # self.CFL = 0.9 / 2.0
        self.dtfixed0 = 1200.0 / self.t_ref
        self.dtfixed = 1200.0 / self.t_ref

        self.inx = 64+1
        self.iny = 1+1
        self.inz = 64+1

        self.recovery_order = RecoveryOrder.SECOND
        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.kp = 0.0
        self.kz = 0.0
        self.km = 0.0
        self.kY = 0.0
        self.kZ = 0.0

        self.tol = 1.e-12
        self.max_iterations = 6000

        self.perturb_type = 'pos_perturb'
        self.blending_mean = 'rhoY' # 1.0, rhoY
        self.blending_conv = 'swe' #theta, rho

        self.continuous_blending = False
        self.no_of_pi_initial = 1
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 0
        self.no_of_hy_transition = 0

        self.blending_weight = 16./16

        self.initial_projection = True
        self.initial_impl_Euler = False

        self.initial_blending = True

        self.column_preconditionr = False
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 0.0

        # self.tout = np.arange(0,2E5+stepsize,stepsize)
        # self.tout = np.arange(0,1E6+10000,10000)[1:]
        stepsize = 86400.0/4
        stepsize = 1200.0 / self.t_ref
        end = 86400.0*3.0 / self.t_ref
        self.tout = np.arange(0,end+stepsize, stepsize)[1:]
        self.stepmax = 2E6

        self.output_base_name = "_swe_dvortex"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.tout[-1])
        
        # aux = 'dvortex_comp_geostrophic'
        aux = 'comp_debug'
        aux = 'psinc_debug'
        aux = 'comp_debug_ib'
        # aux += '_' + self.blending_conv + '_conv'
        # aux += '_' + self.blending_mean + '_mean'
        # aux = 'cb1_w=-6_debug'
        # self.output_suffix += '_w=%i-%i' %(self.blending_weight*16.0,16.0-(self.blending_weight*16.0))
        # aux = 'psinc_bal_debug'
        self.output_suffix = "_%i_%i_%i_%.1f_%s" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1],aux)
        self.aux = aux

        self.stratification = self.stratification_function
        self.rhoe = self.rhoe_function

    def stratification_function(self, y):
        return 1.0

    def rhoe_function(self,rho,u,v,w,p,ud,th):
        Msq = ud.compressibility * ud.Msq
        gm1inv = th.gm1inv

        return p * gm1inv + 0.5 * Msq * rho * (u**2 + v**2 + w**2)

def sol_init(Sol, mpv, elem, node, th, ud, seed=None):
    # wind speed
    u0 = 0.0
    v0 = 0.0
    w0 = 0.0

    igs = elem.igs
    igy = igs[1]

    f0 = ud.coriolis_strength[0]

    i2 = (slice(igs[0],-igs[0]),slice(igs[1],-igs[1]),slice(igs[2],-igs[2]))

    hydrostatic_state(mpv, elem, node, th, ud)

    x, z = elem.x[igs[0]:-igs[0]], elem.z[igs[2]:-igs[2]]
    X, Z = np.meshgrid(x,z)

    Lx = (elem.x[-1] - elem.x[0]) #/ ud.h_ref
    Lz = (elem.z[-1] - elem.z[0]) #/ ud.h_ref
    Lx = 5000.0*1E3 / ud.h_ref #elem.x[-1] - elem.x[0]
    Lz = 4330.0*1E3 / ud.h_ref #elem.z[-1] - elem.z[0]

    Hp = 75.0 / ud.h_ref
    H00 = 450.0 # semi-geostrophic
    H00 = 750.0 # quasi-geostrophic
    # H00 = 10000.0 # incompressible
    if seed is not None:
        H0 = H00 / ud.h_ref + 100.0 * 2.0 * (np.random.random() - 0.5) / ud.h_ref
    else:
        H0 = H00 / ud.h_ref

    g0 = 9.81 / (ud.u_ref / ud.t_ref)
    ud.g0 = g0

    sigx, sigz = sigma(Lx), sigma(Lz)
    
    xc1, xc2 = cs(-1.0, Lx), cs(+1.0, Lx)
    zc1, zc2 = cs(-1.0, Lz), cs(+1.0, Lz)
    xp1, xp2 = ps(X, xc1, Lx, sigx), ps(X, xc2, Lx, sigx)
    zp1, zp2 = ps(Z, zc1, Lz, sigz), ps(Z, zc2, Lz, sigz)
    xpp1, xpp2 = pps(X, xc1, Lx, sigx), pps(X, xc2, Lx, sigx)
    zpp1, zpp2 = pps(Z, zc1, Lz, sigz), pps(Z, zc2, Lz, sigz)

    rho = np.zeros_like((X))
    # initial velocities
    u, v, w = np.zeros_like((X)), 0.0, np.zeros_like((Z))

    exp1 = np.exp(-0.5 * (xp1**2 + zp1**2))
    exp2 = np.exp(-0.5 * (xp2**2 + zp2**2))

    rho[...] = H0 - Hp * ( exp1 + exp2 - (4.0 * np.pi * (sigx * sigz) / (Lx * Lz)) )

    aa = 2
    rho = np.expand_dims(rho, 1)
    rho = np.repeat(rho, elem.icy-aa*igs[1], axis=1)

    u[...] = - g0 * Hp / (f0 * sigz) * (zpp1 * exp1 + zpp2 * exp2).T
    w[...] = + g0 * Hp / (f0 * sigx) * (xpp1 * exp1 + xpp2 * exp2).T 

    Frsq = (u.max()**2 + w.max()**2) / (g0 * rho.max())
    Fr = np.sqrt(Frsq)
    print(Frsq, Fr)
    u = np.expand_dims(u, 1)
    u = np.repeat(u, elem.icy-aa*igs[1], axis=1)
    w = np.expand_dims(w, 1)
    w = np.repeat(w, elem.icy-aa*igs[1], axis=1)

    if (ud.is_compressible) :
        p = g0 / 2.0 * rho**2

        Sol.rho[i2] = rho
        Sol.rhou[i2] = rho * u
        Sol.rhov[i2] = rho * v
        Sol.rhow[i2] = rho * w
        Sol.rhoY[i2] = p**th.gamminv
    else:
        H0 = rho.min()
        min_val = rho.min()

        rho_diff = rho.max() - min_val

        H1 = (rho - min_val)

        p = g0 / 2.0 * H0**2
        H1 = g0 / 2.0 * (H1)**2 

        Sol.rho[i2] = H0
        Sol.rhou[i2] = H0 * u
        Sol.rhov[i2] = H0 * v
        Sol.rhow[i2] = H0 * w
        Sol.rhoY[i2] = H1**th.gamminv
    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    kernel = np.ones((2,2))
    kernel /= kernel.sum()
    if (ud.is_compressible):
        pn = signal.convolve(Sol.rhoY[:,igy,:], kernel, mode='valid')
    else:
        pn = signal.convolve((Sol.rhoY[:,igy,:]), kernel, mode='valid')
        pn -= pn.mean()
        Sol.rhoY[i2] = p**th.gamminv
        set_explicit_boundary_data(Sol,elem,ud,th,mpv)
    
    pn = np.expand_dims(pn, 1)
    pn = np.repeat(pn, node.icy, axis=1)

    mpv.p2_nodes[1:-1,:,1:-1] = pn
    set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    pn = np.expand_dims(mpv.p2_nodes[:,igy,:], 1)
    mpv.p2_nodes[...] = np.repeat(pn[...], node.icy, axis=1)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)
    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)


oo = 0.1

def cs(sgn,L):
    return (0.5 + sgn * oo) * L

def sigma(L):
    return(3.0 / 40 * L)

def ps(xx, xxc, L, sig):
    return L / (np.pi * sig) * np.sin( np.pi / L * (xx - xxc) )

def pps(xx, xxc, L, sig):
    return L / (2.0 * np.pi * sig) * np.sin( 2.0 * np.pi / L * (xx - xxc) )

