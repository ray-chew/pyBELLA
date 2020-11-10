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
    omega = 6.147 * 1E-5  # [rad/s] *
    # omega = 0.7071 * 4. * np.pi
    g0 = 9.81

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
    d_ref = h_ref          # [m]
    t_ref = 1200.0      # [day] -> [s]
    u_ref = h_ref / t_ref
    T_ref = 1.0
    p_ref = 1.0

    R_vap = 1.0
    Q_vap = 1.0

    # u_ref = h_ref / t_ref
    rho_ref = p_ref / (R_gas * T_ref)

    Nsq_ref = 0.0

    i_gravity = np.zeros((3))
    i_coriolis = np.zeros((3))

    tout = np.zeros((2))

    def __init__(self):
        self.h_ref = self.h_ref
        self.d_ref = self.d_ref
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
        self.Msq = 2.0 * self.u_ref * self.u_ref / (self.d_ref * self.g0)
        self.g0 = self.g0 / (self.d_ref / self.t_ref**2)


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
        self.xmax =   5000.0 * 1E3 / self.h_ref
        self.ymin =  -0.5 * self.d_ref
        self.ymax =   0.5 * self.d_ref
        self.zmin =   0.0
        self.zmax =   4330.0 * 1E3 / self.h_ref

        self.u_wind_speed = 0.0
        self.v_wind_speed = 0.0
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

        self.initial_blending = False

        self.initial_projection = True
        self.initial_impl_Euler = False

        self.column_preconditionr = False
        self.synchronize_nodal_pressure = False
        self.synchronize_weight = 0.0

        stepsize = 100
        # self.tout = np.arange(0,2E5+stepsize,stepsize)
        # self.tout = np.arange(0,1E6+100,100)[1:]
        self.tout = np.arange(0,(86400.0*10.0+1200.0) / self.t_ref,1200.0 / self.t_ref)[1:]
        self.stepmax = 2E6

        self.output_base_name = "_swe"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.tout[-1])
        
        # aux = 'icshear_pi'
        aux = 'icshear_comp'
        aux = 'debug_psinc'
        # aux += '_' + self.blending_conv + '_conv'
        # aux += '_' + self.blending_mean + '_mean'
        # aux = 'cb1_w=-6_debug'
        # self.output_suffix += '_w=%i-%i' %(self.blending_weight*16.0,16.0-(self.
    # kappa = 0.05
    # rho = 256. / N * Lx*Lz / g0 * ((1.0 / np.pi * (Z - np.pi) * np.exp(-2.0 * (Z - np.pi)**2) * (1.0 + kappa * np.sin(2.0 * X)) + 1.0 )**(-1) - 1.0) + 1.0blending_weight*16.0))
        # aux = 'psinc_bal_debug'
        self.output_suffix = "_%i_%i_%i_%.1f_%s" %(self.inx-1,self.iny-1,self.inz-1,self.tout[-1],aux)

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

    # initial velocities
    u, v, w = 0.0, 0.0 , 0.0

    igs = elem.igs
    igy = igs[1]

    g0 = ud.g0

    f0 = ud.coriolis_strength[0]
    # ud.Msq = 1e-16

    i2 = (slice(igs[0],-igs[0]),slice(igs[1],-igs[1]),slice(igs[2],-igs[2]))
    i2 = (slice(None,),slice(None,),slice(None,))

    # x, z = elem.x[igs[0]:-igs[0]], elem.z[igs[2]:-igs[2]]
    x, z = elem.x, elem.z
    X, Z = np.meshgrid(x,z)

    # hydrostatic_state(mpv, elem, node, th, ud)

    # Lx = ud.xmax
    # Lz = ud.zmax
    Lx = elem.x[-1] - elem.x[0]
    Lz = elem.z[-1] - elem.z[0]
    Lx = 5000.0*1E3 / ud.h_ref #elem.x[-1] - elem.x[0]
    Lz = 4330.0*1E3 / ud.h_ref #elem.z[-1] - elem.z[0]

    sigz = 1.0 / 12
    lambdx = 0.5
    kappa = 0.1

    xp = X / Lx

    Hp = 30.0 / ud.d_ref

    if seed is not None:
        H0 = 1076.0 / ud.d_ref + 10.0 * np.random.random() / ud.h_ref
        # print(H0)
    else:
        H0 = 1076.0 / ud.d_ref

    exp0 = np.exp(- zp(Z, Lz)**2 / (2.0 * sigz**2) + 0.5)

    rho = H0 - Hp * zpp(Z, Lz) / sigz * exp0 * (1.0 + kappa * np.sin(2.0 * np.pi * xp / lambdx))

    u = g0 * Hp / (sigz * f0 * Lz) * (cz(Z,Lz) - zpp(Z, Lz)**2 / sigz**2) * exp0 * (1.0 + kappa * np.sin(2.0 * np.pi * xp / lambdx))

    w = -g0 * Hp / (f0 * Lx) * 2.0 * np.pi * kappa * zpp(Z,Lz) / (lambdx * sigz) * exp0 * np.cos(2.0 * np.pi * xp / lambdx)

    # rho, u, w = rho.T, u.T[::-1,:], w.T[:,::-1]
    rho, u, w = rho.T, u.T, w.T

    aa = 2

    rho = np.expand_dims(rho, 1)
    rho = np.repeat(rho, elem.icy-aa*igs[1], axis=1)
    u = np.expand_dims(u, 1)
    u = np.repeat(u, elem.icy-aa*igs[1], axis=1)
    w = np.expand_dims(w, 1)
    w = np.repeat(w, elem.icy-aa*igs[1], axis=1)

    fr = ((u.max()*ud.u_ref)**2 + (w.max()*ud.u_ref)**2)**0.5 / (g0 * (ud.u_ref / ud.t_ref) * (rho.max() * (ud.h_ref)))**0.5
    print("Froude number =", fr)
    print("Fr2 = ", fr**2)

    if (ud.is_compressible):
        Sol.rho[i2] = rho
        Sol.rhou[i2] = rho * u
        Sol.rhov[i2] = rho * v
        Sol.rhow[i2] = rho * w
        Sol.rhoY[i2] = rho
    else:
        H0 = rho.mean()
        min_val = rho.mean()

        H1 = (rho - min_val)

        Sol.rho[i2] = H0
        Sol.rhou[i2] = H0 * u
        Sol.rhov[i2] = H0 * v
        Sol.rhow[i2] = H0 * w
        Sol.rhoY[i2] = H1

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    kernel = np.ones((2,2))
    kernel /= kernel.sum()

    if (ud.is_compressible):
        pn = (Sol.rho[:,igy,:] - H0)
        pn /= ud.Msq
        pn = signal.convolve( pn, kernel, mode='valid')
    else:
        pn = Sol.rhoY[:,igy,:] #/ ud.Msq
        pn = signal.convolve(pn, kernel, mode='valid')
        Sol.rhoY[i2] = H0
        set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    pn = np.expand_dims(pn, 1)
    pn = np.repeat(pn, node.icy, axis=1)

    # mpv.p2_nodes[...] = pn
    mpv.p2_nodes[1:-1,:,1:-1] = pn
    set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    hydrostatic_state(mpv, elem, node, th, ud)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)

def icshear(X,Y,f0):
    Y = Y - np.pi
    dq = 1./np.pi * Y * np.exp(-2. * Y**2) * (1. + 0.1 * np.sin(2. * X))
    
    dq -= np.mean(dq)
    
    q = f0 * (1. + dq)
    
    H = f0/q
    H = 1. + (H - np.mean(H))
    return H

def interpolate(x,y,Un,N,L,Xn=None,Yn=None):
    Nint = int(N)
    x_extended = np.arange(Nint+1)/N*L
    y_extended = np.copy(x_extended)
    
    Un_top = np.hstack((Un[:,:],Un[:,0].reshape(-1,1)))
    Un_bottom = np.hstack((Un[0,:],[Un[0,0]]))
    Un = np.vstack((Un_top,Un_bottom))
    
    if ((Xn is None) and (Yn is None)):
        U = RectBivariateSpline(y_extended,x_extended,Un)
    else:
        U = RectBivariateSpline(y_extended,x_extended,Un)(Yn.ravel(),Xn.ravel(),grid=False).reshape(Nint,Nint)
        
    return U

def Dx(u,N,L,I):
    us = np.roll(u,-I,axis=1)
    return I * N/L * (us-u)

def Dy(u,N,L,I):
    us = np.roll(u,-I,axis=0)
    return I * N/L * (us-u)



def xp(x, Lx):
    return x / Lx

def cz(z, Lz):
    cc = z - (Lz / 2)
    return np.cos( 2.0 * np.pi / Lz * cc )

def zp(z, Lz):
    cc = z - (Lz / 2)
    return ( 1.0 / np.pi * np.sin(np.pi / Lz * cc) )

def zpp(z, Lz):
    cc = z - (Lz / 2)
    return ( 1.0 / (2.0 * np.pi) * np.sin(2.0 * np.pi / Lz * cc) )