import numpy as np
from scipy import signal
from inputs import boundary

class Blend(object):
    def __init__(self,ud):
        self.bb = False
        self.cb = ud.continuous_blending
        self.psinc_init = ud.no_of_pi_initial
        self.psinc_trans = ud.no_of_pi_transition

        if self.psinc_init > 0 and self.psinc_trans == 0:
            self.bb = True

        self.c_init = self.criterion_init
        self.c_trans = self.criterion_trans

        self.fac = ud.Msq
    
    def criterion_init(self, step):
        return step == (self.psinc_init) and self.cb and self.bb

    def criterion_trans(self, step):
        return step <= self.psinc_trans and self.cb and not self.bb
        
    def convert_p2n(self, p2n):
        ndim =  p2n.ndim
        dp2n = p2n - p2n.mean()

        self.kernel = np.ones([2]*ndim)
        dp2c = signal.fftconvolve(dp2n,self.kernel, mode='valid') / self.kernel.sum()

        self.dp2n = dp2n - dp2n.mean()
        self.dp2c = dp2c - dp2c.mean()

        return dp2c

    def update_Sol(self, Sol, elem, node, th, ud, mpv, label=None,writer=None):
        # if writer != None: writer.populate(str(label)+'_before_blending', 'dp2n', self.dp2n)
        # if writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_before_blending')

        rho = np.copy(Sol.rho)
        rhoY = np.copy(Sol.rhoY)

        Y = rhoY / rho

        # rhoYc = (Sol.rhoY**th.gm1 + self.fac * self.dp2c)**(th.gm1inv)
        
        rhoYc = (1.0 + self.fac * self.dp2c)**(th.gm1inv)

        alpha = rhoYc / Sol.rhoY

        ### keep theta, convert rho

        Sol.rho[...] = rho*alpha
        Sol.rhoY[...] = (Sol.rho * Y)

        rho_fac = Sol.rho / rho
        Sol.rhou[...] *= rho_fac
        Sol.rhov[...] *= rho_fac
        Sol.rhow[...] *= rho_fac
        Sol.rhoX[...] *= rho_fac

        ### keep rho, convert theta

        # Yc = Y * alpha
        # Sol.rhoY[...] = rho * Yc
        # Sol.rhoX[...] = rho / Yc

        # if writer != None: writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_blending')

    def update_p2n(self,Sol, mpv, node, th, ud):
        # rhoYncomp = np.zeros_like(self.dp2n)
        # rhoYncomp[1:-1,1:-1] = signal.fftconvolve(Sol.rhoY,self.kernel,mode='valid') / self.kernel.sum()
        # boundary.set_ghostnodes_p2(rhoYncomp,node,ud)
        # mpv.p2_nodes[...] = (rhoYncomp**th.gm1) - 1.0 + self.dp2n
        # mpv.p2_nodes -= mpv.p2_nodes.mean()

        mpv.p2_nodes = self.dp2n

        ### shift p2n
        # tol = .001
        # diff = mpv.p2_nodes - 0.0
        # locs = np.where(np.abs(diff) < tol)
        # mpv.p2_nodes[locs] -= diff[locs]