import numpy as np
import dycore.physics.hydrostatics as hydrostatic
import dycore.utils.boundary as bdry

class UserData(object):
    # Nsq_ref = grav * 1.3e-05

    def __init__(self):
        self.grav = 10.0             # [m/s^2]
        self.t_ref = 1000.0           # [s]

        ##########################################
        # NUMERICS
        ##########################################
        self.CFL  = 0.5
        self.dtfixed0 = 100.0
        self.dtfixed = 100.0

        self.inx = 160+1
        self.iny = 80+1
        self.inz = 1

        

        self.tout = np.arange(0.0,1.01,0.01)[10:]
        self.stepmax = 10000

        self.output_base_name = "_rising_bubble"
        # if self.is_compressible == 1:
        #     self.output_suffix = "_%i_%i_%.1f_comp" %(self.inx-1,self.iny-1,self.tout[-1])
        # if self.is_compressible == 0:
        #     self.output_suffix = "_%i_%i_%.1f_psinc" %(self.inx-1,self.iny-1,self.tout[-1])
        # if self.continuous_blending == True:
        #     self.output_suffix = "_%i_%i_%.1f" %(self.inx-1,self.iny-1,self.tout[-1])
        
        aux = 'debug_imbal_CFLfixed'
        self.aux = aux
        # self.output_suffix = "_%i_%i_%.1f_%s" %(self.inx-1,self.iny-1,self.tout[-1],aux)
        # self.output_suffix += '_w=%i-%i' %(self.blending_weight*16.0,16.0-(self.blending_weight*16.0))


def sol_init(Sol, mpv, elem, node, th, ud, seed=None):
    u0 = ud.u_wind_speed
    v0 = ud.v_wind_speed
    w0 = ud.w_wind_speed
    delth = 2.0         # [K]
    
    y0 = 0.2
    r0 = 0.2

    g = ud.gravity_strength[1]
    # print(ud.rho_ref)

    hydrostatic.state(mpv, elem, node, th, ud)

    x = elem.x
    y = elem.y

    x, y = np.meshgrid(x,y)

    if seed != None:
        np.random.seed(seed)
        # y0 += (np.random.random()-.5)/2.0
        # delth += 10.0*(np.random.random()-.5)
        delth += 10.0*(np.random.random())
    
    if 'truth' in ud.aux:
        np.random.seed(1234)
        # delth += 10.0*(np.random.random()-.5)
        delth += 10.0*(np.random.random())
    print(delth)
    
    r = np.sqrt((x)**2 + (y-y0)**2) / r0

    p = np.repeat(mpv.HydroState.p0.reshape(1,-1),elem.icx,axis=0)
    rhoY = np.repeat(mpv.HydroState.rhoY0.reshape(1,-1),elem.icx,axis=0)

    perturbation = (delth/300.0) * (np.cos(0.5 * np.pi * r)**2)
    perturbation[np.where(r > 1.0)] = 0.0
    rho = rhoY / (ud.stratification(y) + perturbation.T)

    x_idx = slice(None)
    y_idx = slice(None)

    u, v, w = u0, v0, w0

    Sol.rho[x_idx,y_idx] = rho
    Sol.rhou[x_idx,y_idx] = rho * u
    Sol.rhov[x_idx,y_idx] = rho * v
    Sol.rhow[x_idx,y_idx] = rho * w
    Sol.rhoY[x_idx,y_idx] = rhoY

    p = mpv.HydroState_n.p0[0]
    rhoY = mpv.HydroState_n.rhoY0[0]
    mpv.p2_nodes[...] = (p - mpv.HydroState_n.p0[0]) / rhoY / ud.Msq

    bdry.set_explicit_boundary_data(Sol,elem,ud,th,mpv)

    return Sol