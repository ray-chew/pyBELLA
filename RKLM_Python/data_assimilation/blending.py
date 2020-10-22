import numpy as np
from scipy import signal
from inputs import boundary

from termcolor import colored
from copy import deepcopy
class Blend(object):
    """
    Class that takes care of the blending interface.
    """

    def __init__(self,ud):
        self.bb = False
        self.cb = ud.continuous_blending
        self.psinc_init = ud.no_of_pi_initial
        self.psinc_trans = ud.no_of_pi_transition

        if self.psinc_init > 0 and self.psinc_trans == 0 and self.cb:
            self.bb = True

        self.c_init = self.criterion_init
        self.c_trans = self.criterion_trans

        self.fac = ud.Msq
    
    def criterion_init(self, step):
        return step == (self.psinc_init) and self.cb and self.bb

    def criterion_trans(self, step):
        return step <= self.psinc_trans and self.cb and not self.bb
        
    def convert_p2n(self, p2n):
        ndim = p2n.ndim
        dp2n = p2n - p2n.mean()

        self.kernel = np.ones([2]*ndim)
        dp2c = signal.fftconvolve(dp2n,self.kernel, mode='valid') / self.kernel.sum()

        self.dp2n = dp2n - dp2n.mean()
        self.dp2c = dp2c - dp2c.mean()

        return dp2c

    def update_Sol(self, Sol, elem, node, th, ud, mpv, sgn, label=None,writer=None):
        if writer != None: writer.populate(str(label)+'_before_blending', 'dp2n', self.dp2n)
        if writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_blending')

        if sgn == 'bef':
            sign = -1.0
        elif sgn == 'aft':
            sign = +1.0
        else:
            assert 0, "sgn == bef or sgn == aft"

        rho = np.copy(Sol.rho)
        rhoY = np.copy(Sol.rhoY)

        Y = rhoY / rho

        if ud.blending_mean == 'rhoY':
            rhoYc = (rhoY**th.gm1 + sign * self.fac * self.dp2c)**(th.gm1inv)
        elif ud.blending_mean == '1.0':
            rhoYc = (1.0 + sign * self.fac * self.dp2c)**(th.gm1inv)

        alpha = rhoYc / Sol.rhoY

        if ud.blending_conv == 'rho':
        ### keep theta, convert rho
            Sol.rho[...] = rho*alpha
            Sol.rhoY[...] = (Sol.rho * Y)

            rho_fac = Sol.rho / rho
            Sol.rhou[...] *= rho_fac
            Sol.rhov[...] *= rho_fac
            Sol.rhow[...] *= rho_fac
            Sol.rhoX[...] *= rho_fac

        elif ud.blending_conv == 'theta':
        ### keep rho, convert theta
            Yc = Y * alpha
            Sol.rhoY[...] = rho * Yc
            Sol.rhoX[...] = rho / Yc
        else:
            assert(0, "ud.blending_conv undefined.")

        if writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_blending')

    def update_p2n(self,Sol, mpv, node, th, ud):
        mpv.p2_nodes = self.dp2n


def do_comp_to_psinc_conv(Sol, mpv, bld, elem, node, th, ud, label, writer):
    print(colored("Converting COMP to PSINC",'blue'))
    dp2n = mpv.p2_nodes
    bld.convert_p2n(dp2n)
    bld.update_Sol(Sol,elem,node,th,ud,mpv,'bef',label=label,writer=writer)
    bld.update_p2n(Sol,mpv,node,th,ud)


def do_psinc_to_comp_conv(Sol, flux, mpv, bld, elem, node, th, ud, label, writer, time_update, step, window_step, t, dt):
    print(colored("Blending... step = %i" %window_step,'blue'))
    Sol_freeze = deepcopy(Sol)
    mpv_freeze = deepcopy(mpv)

    ret = time_update(Sol,flux,mpv, t, t+dt, ud, elem, node, [0,step-1], th, bld=None, writer=None, debug=False)

    fac_old = ud.blending_weight
    fac_new = 1.0 - fac_old
    dp2n = (fac_new * ret[2].p2_nodes + fac_old * mpv_freeze.p2_nodes)

    if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_start', mpv_freeze.p2_nodes)
    if writer != None: writer.populate(str(label)+'_after_full_step', 'p2_end', ret[2].p2_nodes)
    Sol = Sol_freeze
    mpv = mpv_freeze

    if writer != None: writer.populate(str(label)+'_after_full_step', 'dp2n', dp2n)
    print(colored("Converting PSINC to COMP",'blue'))
    bld.convert_p2n(dp2n)
    bld.update_Sol(Sol,elem,node,th,ud,mpv,'aft',label=label,writer=writer)
    bld.update_p2n(Sol,mpv,node,th,ud)


def do_swe_to_lake_conv(Sol, mpv, elem, node, ud, th, writer, label, debug):
    print(colored("swe to lake conversion...",'blue'))

    H1 = deepcopy(Sol.rho[:,2,:])
    setattr(ud,'mean_val',H1.mean())
    H1 = (H1 - ud.mean_val)

    # pn = mpv.p2_nodes[1:-1,2,1:-1]

    # kernel = np.ones((2,2))
    # kernel /= kernel.sum()

    # do node-to-cell averaging
    # pn = signal.convolve(H1, kernel, mode='valid')

    Sol.rhou[...] = Sol.rhou / Sol.rho * ud.mean_val
    Sol.rhov[...] = Sol.rhov / Sol.rho * ud.mean_val
    Sol.rhow[...] = Sol.rhow / Sol.rho * ud.mean_val
    Sol.rhoY[...] = Sol.rhoY / Sol.rho * ud.mean_val
    Sol.rho[...] = ud.mean_val

    # pn = np.expand_dims(pn, axis=1)
    # mpv.p2_nodes[1:-1,:,1:-1] = np.repeat(pn[...], node.icy, axis=1)
    # boundary.set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_swe_to_lake')


def do_lake_to_swe_conv(Sol, mpv, elem, node, ud, th, writer, label, debug):
    if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_lake_time_step')
    print(colored("lake to swe conversion...",'blue'))

    H10 = deepcopy(mpv.p2_nodes[:,2,:])
    H10 -= H10.mean()

    # define 2D kernel
    kernel = np.ones((2,2))
    kernel /= kernel.sum()

    # do node-to-cell averaging
    H1 = signal.convolve(H10, kernel, mode='valid')
    H1 = ud.mean_val + ud.Msq*H1
    print(colored(H1.max(), 'red'))

    # project H1 back to horizontal slice with ghost cells
    H1 = np.expand_dims(H1, axis=1)
    H1 = np.repeat(H1, elem.icy, axis=1)

    Sol.rho[...] = H1
    Sol.rhou[...] = Sol.rhou / ud.mean_val * Sol.rho
    Sol.rhov[...] = Sol.rhov / ud.mean_val * Sol.rho
    Sol.rhow[...] = Sol.rhow / ud.mean_val * Sol.rho
    Sol.rhoY[...] = Sol.rhoY / ud.mean_val * Sol.rho
    
    pn = H10[1:-1,1:-1]
    pn = np.expand_dims(pn, axis=1)
    mpv.p2_nodes[1:-1,:,1:-1] = np.repeat(pn[...], node.icy, axis=1)
    boundary.set_ghostnodes_p2(mpv.p2_nodes,node,ud)

    if debug == True: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_lake_to_swe')


