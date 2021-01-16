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
    omega = 0.7071 * 4. * np.pi 

    R_gas = 1.0
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
    T_ref = 1.0
    p_ref = 1e5

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

        self.g = 9.81 * self.h_ref / (self.R_gas * self.T_ref)
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

        # self.xmin = - 5E+5
        # self.xmax =   5E+5
        # self.ymin = - 5E+5
        # self.ymax =   5E+5
        self.xmin = - 0.0
        self.xmax =   2.0 * np.pi
        self.ymin = - 0.0
        self.ymax =   2.0 * np.pi
        self.zmin = - 0.5
        self.zmax =   0.5

        self.wind_speed = 0.0

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
        self.dtfixed0 = 1.0/24.*3
        self.dtfixed = 1.0/24.*3

        self.inx = 128+1
        self.iny = 128+1
        self.inz = 1

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

        stepsize = 100
        # self.tout = np.arange(0,2E5+stepsize,stepsize)
        self.tout = np.arange(0,1E6+100,100)#[1:]
        self.stepmax = 201

        self.output_base_name = "_swe"
        if self.is_compressible == 1:
            self.output_suffix = "_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.is_compressible == 0:
            self.output_suffix = "_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        if self.continuous_blending == True:
            self.output_suffix = "_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.tout[-1])
        
        aux = 'icshear'
        # aux += '_' + self.blending_conv + '_conv'
        # aux += '_' + self.blending_mean + '_mean'
        # aux = 'cb1_w=-6_debug'
        # self.output_suffix += '_w=%i-%i' %(self.blending_weight*16.0,16.0-(self.blending_weight*16.0))
        # aux = 'psinc_bal_debug'
        self.output_suffix = "_%i_%i_%.1f_%s" %(self.inx-1,self.iny-1,self.tout[-1],aux)

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

    g0 = ud.gravity_strength[1]
    f0 = ud.coriolis_strength[0]
    ud.Msq = 1.0

    i2 = (slice(igs[0],-igs[0]),slice(igs[1],-igs[1]))

    hydrostatic_state(mpv, elem, node, th, ud)

    X,Y = np.meshgrid(elem.x,elem.y)

    Lx = elem.x[-1] - elem.x[0]
    Ly = elem.y[-1] - elem.y[0]
    N = len(elem.x) ## assuming square grid
    dx = np.diff(elem.x)[0] ## assuming dx == dy

    factor = 256./N * np.sqrt(3.)
    g0 = factor**2 * Lx*Ly
    # f0 = 0.7071 * 4. * np.pi ## temp
    rho = icshear(X,Y,f0)

    rho = 256. / N * Lx*Ly / g0 * (rho-np.mean(rho)) + np.mean(rho)

    Hx = Dx(rho,N,Lx,-1)
    Hy = Dy(rho,N,Ly,-1)

    Uv = -g0 * Hy / f0
    Vu = +g0 * Hx / f0

    u = interpolate(elem.x,elem.y,Uv,N,Lx,X,Y+dx/2)
    v = interpolate(elem.x,elem.y,Vu,N,Ly,X+dx/2.,Y)

    rho[i2] = np.load('H0.npy')
    u[i2] = np.load('Uc0.npy')
    v[i2] = np.load('Vc0.npy')    

    p = g0 / 2.0 * rho**2

    Sol.rho[...] = rho
    Sol.rhou[...] = rho * u
    Sol.rhov[...] = rho * v
    Sol.rhow[...] = rho * w

    if (ud.is_compressible) :
        Sol.rhoY[...] = p**th.gamminv
    else:
        Sol.rhoY[...] = 1.0
    set_explicit_boundary_data(Sol,elem,ud,th,mpv)
    
    # rho = np.pad(rho,2,mode='wrap')

    # kernel = np.ones((2,2))
    # kernel /= kernel.sum()
    # rho_n = signal.convolve(Sol.rho, kernel, mode='valid')

    points = np.zeros((Sol.rhoY[...].flatten().shape[0],2))
    points[:,0] = X[...].flatten()
    points[:,1] = Y[...].flatten()

    values = (p**th.gamminv).flatten()

    grid_x, grid_y = np.meshgrid(node.x,node.y)

    from scipy.interpolate import griddata
    pn = griddata(points, values, (grid_x, grid_y), method='cubic')

    mpv.p2_nodes[...] = pn
    set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    ud.nonhydrostasy = float(ud.is_nonhydrostatic)
    ud.compressibility = float(ud.is_compressible)

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    if ud.initial_projection == False:
        is_compressible = np.copy(ud.is_compressible)
        compressibility = np.copy(ud.compressibility)
        ud.is_compressible = 0
        ud.compressibility = 0.0

        p2aux = np.copy(mpv.p2_nodes)

        Sol.rhou -= u0 * Sol.rho
        Sol.rhov -= v0 * Sol.rho

        euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, 0.0, ud.dtfixed, 0.5)
        # euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, 0.0, 2.5e-3, 0.5)

        mpv.p2_nodes[...] = p2aux
        mpv.dp2_nodes[...] = 0.0

        Sol.rhou += u0 * Sol.rho
        Sol.rhov += v0 * Sol.rho

        ud.is_compressible = is_compressible
        ud.compressibility = compressibility

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