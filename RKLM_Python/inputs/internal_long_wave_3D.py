import numpy as np
from inputs.enum_bdry import BdryType
from management.enumerator import TimeIntegrator, MolecularTransport,HillShapes, BottomBC, LimiterType, RecoveryOrder
from physics.hydrostatics import hydrostatic_state, hydrostatic_column, hydrostatic_initial_pressure
from inputs.boundary import set_explicit_boundary_data, set_ghostcells_p2, set_ghostnodes_p2
from management.variable import Vars, States
from physics.low_mach.mpv import MPV
from discretization import kgrid

class UserData(object):
    NSPEC = 1
    BOUY = 0

    grav = 9.81
    omega = -50.0*0.0001

    R_gas = 287.4
    R_vap = 461.0
    Q_vap = 2.53e+06
    gamma = 1.4

    h_ref = 10000.0
    t_ref = 100.0
    T_ref = 300.00
    p_ref = 1e+5
    u_ref = h_ref / t_ref
    rho_ref = p_ref / (R_gas * T_ref)

    Nsq_ref = 1.0e-4

    # planetary -> 160.0;  long-wave -> 20.0;  standard -> 1.0;
    scale_factor = 20.0

    i_gravity = np.zeros((3))
    i_coriolis = np.zeros((3))

    tout = np.zeros((2))

    def __init__(self):
        self.scale_factor = self.scale_factor

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
        self.is_compressible = 0
        self.is_ArakawaKonor = 0

        self.compressibility = 0.0
        self.acoustic_timestep = 0
        self.Msq = self.u_ref * self.u_ref / (self.R_gas * self.T_ref)

        self.gravity_strength = np.zeros((3))
        self.coriolis_strength = np.zeros((3))

        self.gravity_strength[1] = self.grav * self.h_ref / (self.R_gas * self.T_ref)
        # self.coriolis_strength[0] = self.omega * self.t_ref
        self.coriolis_strength[2] = self.omega * self.t_ref
        # self.coriolis_strength[1] = self.omega * self.t_ref

        for i in range(3):
            if (self.gravity_strength[i] > np.finfo(np.float).eps) or (i == 1):
                self.i_gravity[i] = 1
                self.gravity_direction = i

            if (self.coriolis_strength[i] > np.finfo(np.float).eps):
                self.i_coriolis[i] = 1

        self.xmin = -15.0 * self.scale_factor
        self.xmax =  15.0 * self.scale_factor
        self.ymin =   0.0
        self.ymax =   1.0
        self.zmin = - 1.0
        self.zmax =   1.0

        self.u_wind_speed = 0.0 * 20.0 / self.u_ref
        self.v_wind_speed = 0.0
        self.w_wind_speed = 0.0

        self.bdry_type = np.empty((3), dtype=object)
        self.bdry_type[0] = BdryType.PERIODIC
        self.bdry_type[1] = BdryType.WALL
        self.bdry_type[2] = BdryType.PERIODIC

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL  = 0.9

        self.dtfixed0 = 10.0 * ( 12.5 / 15.0) * 0.5 * self.scale_factor * 30.0 / self.t_ref
        self.dtfixed = 10.0 * (12.5 / 15.0) * 0.5 * self.scale_factor * 30.0 / self.t_ref

        self.dtfixed0 = 1.0
        self.dtfixed = 1.0

        self.inx = 301+1
        self.iny = 10+1
        self.inz = 1+1

        self.limiter_type_scalars = LimiterType.NONE
        self.limiter_type_velocity = LimiterType.NONE

        self.initial_projection = False

        self.tout =  [self.scale_factor * 1.0 * 3000.0 / self.t_ref]
        # self.tout = np.arange(0,self.tout[0]+60.0,60.0)[1:]

        self.tol = 1.e-16
        self.stepmax = 100000
        self.max_iterations = 6000

        self.continuous_blending = False
        self.no_of_pi_initial = 0
        self.no_of_pi_transition = 0
        self.no_of_hy_initial = 1
        self.no_of_hy_transition = 0

        self.blending_weight = 0./16
        self.blending_mean = 'rhoY' # 1.0, rhoY
        self.blending_conv = 'rho' #theta, rho
        self.blending_type = 'half' # half, full

        self.initial_blending = False

        self.output_base_name = "_internal_long_wave"
        if self.scale_factor < 10.0:
            self.scale_tag = "standard"
        elif self.scale_factor == 20.0:
            self.scale_tag = "long"
        elif self.scale_factor == 160.0:
            self.scale_tag = "planetary"
        else:
            assert(0, "scale factor unsupported")

        if self.is_nonhydrostatic and self.is_compressible:
            self.h_tag = 'nonhydro'
        elif self.is_nonhydrostatic and not self.is_compressible:
            self.h_tag = 'psinc'
        elif not self.is_nonhydrostatic:
            self.h_tag = 'hydro'

        aux = 'vcor'
        if len(aux) > 0:
            aux = self.scale_tag + "_" + self.h_tag + "_" + aux
        else:
            aux = self.scale_tag + "_" + self.h_tag
        self.aux = aux

        self.output_timesteps = True

        self.stratification = self.stratification_function
        self.molly = self.molly_function
        self.rhoe = self.rhoe_method
        
    def stratification_function(self, y):
        Nsq = self.Nsq_ref * self.t_ref * self.t_ref
        g = self.gravity_strength[1] / self.Msq

        return np.exp(Nsq * y / g)

    def molly_function(self, x):
        del0 = 0.25
        L = self.xmax - self.xmin
        xi_l = np.minimum(1.0, (x - self.xmin) / (del0 * L))
        xi_r = np.minimum(1.0, (self.xmax - x) / (del0 * L))

        return 0.5 * np.minimum(1.0 - np.cos(np.pi * xi_l), 1.0 - np.cos(np.pi * xi_r))

    @staticmethod
    def rhoe_method(rho,u,v,w,p,ud,th):
        Msq = ud.compressibility * ud.Msq

        gm1inv = th.gm1inv
        return (p * gm1inv + 0.5 * Msq * rho * (u*u + v*v + w*w))

def sol_init(Sol, mpv, elem, node, th, ud, seeds=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed
    delth = 0.01 / ud.T_ref
    xc = -1.0 * ud.scale_factor * 50.0e+3 / ud.h_ref
    xc = 0.0
    a = ud.scale_factor * 5.0e+3 / ud.h_ref

    hydrostatic_state(mpv, elem, node, th, ud)

    # generate sizes of data containers assuming
    # they are true vertical slices.
    g_slc = kgrid.Grid(ud.inx, ud.iny, 1, ud.xmin, ud.xmax, ud.ymin, ud.ymax, 0, 1)
    elem_slc = kgrid.ElemSpaceDiscr(g_slc)
    node_slc = kgrid.NodeSpaceDiscr(g_slc)

    HySt  = States(node_slc.sc,ud)
    HyStn = States(node_slc.sc,ud)

    x = elem.x.reshape(-1,1)
    y = elem.y.reshape(1,-1)
    Y = ud.stratification(y) + delth * ud.molly(x) * np.sin(np.pi * y) / (1.0 + (x-xc)**2 / (a**2))

    xn = node.x[:-1].reshape(-1,1)
    yn = node.y[:-1].reshape(1,-1)
    Yn = ud.stratification(yn) + delth * ud.molly(xn) * np.sin(np.pi * yn) / (1.0 + (xn-xc)**2 / (a**2))

    hydrostatic_column(HySt, HyStn, Y, Yn, elem_slc, node_slc, th, ud)

    x_idx = slice(None)
    y_idx = slice(elem.igy,-elem.igy+1)
    xc_idx = slice(0,-1)
    yc_idx = slice(0,-1)
    c_idx = (xc_idx, yc_idx)

    Sol_slc = Vars(elem_slc.sc, ud)
    mpv_slc = MPV(elem_slc, node_slc, ud)

    hydrostatic_state(mpv_slc, elem_slc, node_slc, th, ud)

    u, v, w = u0, v0, w0
    if ud.is_compressible:
        p = HySt.p0[:,y_idx][c_idx]
        rhoY = HySt.rhoY0[:,y_idx][c_idx]
    else:
        p = mpv_slc.HydroState.p0[y_idx]
        rhoY = mpv_slc.HydroState.rhoY0[y_idx]

    rho = rhoY / Y[:,y_idx]
    Sol_slc.rho[x_idx,y_idx] = rho
    Sol_slc.rhou[x_idx,y_idx] = rho * u
    Sol_slc.rhov[x_idx,y_idx] = rho * v
    Sol_slc.rhow[x_idx,y_idx] = rho * w
    Sol_slc.rhoe[x_idx,y_idx] = ud.rhoe(rho, u ,v, w, p, ud, th)
    Sol_slc.rhoY[x_idx,y_idx] = rhoY

    mpv_slc.p2_cells[x_idx,y_idx] = HySt.p20[x_idx,y_idx][c_idx]

    Sol_slc.rhoX[x_idx,y_idx] = Sol_slc.rho[x_idx,y_idx] * (1.0 / Y[:, y_idx] - mpv.HydroState.S0[y_idx])

    mpv_slc.p2_nodes[:,elem.igy:-elem.igy] = HyStn.p20[:,elem.igy:-elem.igy]

    ud.nonhydrostasy = 1.0 if ud.is_nonhydrostatic == 1 else 0.0
    ud.compressibility = 1.0 if ud.is_compressible == 1 else 0.0

    if 'imbal' in ud.aux:
        mpv.p2_nodes[...] = 0.0
    
    hydrostatic_initial_pressure(Sol_slc,mpv_slc,elem_slc,node_slc,ud,th)

    Sol.rho[...] = Sol_slc.rho[...,np.newaxis]
    Sol.rhou[...] = Sol_slc.rhou[...,np.newaxis]
    Sol.rhov[...] = Sol_slc.rhov[...,np.newaxis]
    Sol.rhoY[...] = Sol_slc.rhoY[...,np.newaxis]
    Sol.rhoX[...] = Sol_slc.rhoX[...,np.newaxis]
    mpv.p2_nodes[...] = mpv_slc.p2_nodes[...,np.newaxis]

    set_explicit_boundary_data(Sol,elem,ud,th,mpv)
    set_ghostnodes_p2(mpv.p2_nodes, node, ud)

    return Sol

def T_from_p_rho(p, rho):
    return np.divide(p,rho)