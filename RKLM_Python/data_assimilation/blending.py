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
        return step == self.psinc_init and self.cb and self.bb

    def criterion_trans(self, step):
        return step <= self.psinc_trans and self.cb and not self.bb

    def interface(self):
        None
        
    def convert_p2n(self, p2n, p2n0):
        self.ps = p2n0 - p2n0.mean()
        self.pn = p2n - p2n.mean()

        ndim =  self.pn.ndim

        # if self.cb and not self.bb:
        #     dp2n = p2n0
        # else:
        #     dp2n = p2n0
        dp2n = 0.5*(self.ps + self.pn)
            

        self.kernel = np.ones([2]*ndim)
        dp2c = signal.fftconvolve(dp2n,self.kernel, mode='valid') / self.kernel.sum()

        self.pnc = signal.fftconvolve(self.pn,self.kernel, mode='valid') / self.kernel.sum()

        self.psc = signal.fftconvolve(self.ps,self.kernel, mode='valid') / self.kernel.sum()

        self.dp2n = dp2n - dp2n.mean()
        self.dp2c = dp2c - dp2c.mean()

        return dp2c

    def update_Sol(self, Sol, th, ud):
        rho = np.copy(Sol.rho)
        rhoY = np.copy(Sol.rhoY)

        Y = rhoY / rho

        Sol.rhoY = (Sol.rhoY**th.gm1 + self.fac * self.dp2c)**(th.gm1inv)

        Sol.rho[...] = Sol.rhoY / Y

        rho_fac = Sol.rho / rho
        Sol.rhou[...] *= rho_fac
        Sol.rhov[...] *= rho_fac
        Sol.rhow[...] *= rho_fac
        Sol.rhoX[...] *= rho_fac

    def update_p2n(self,Sol, mpv, node, th, ud):
        # rhoYncomp = np.zeros_like(self.dp2n)
        # rhoYncomp[1:-1,1:-1] = signal.fftconvolve(Sol.rhoY,self.kernel,mode='valid') / self.kernel.sum()
        # boundary.set_ghostnodes_p2(rhoYncomp,node,ud)
        # mpv.p2_nodes[...] = (rhoYncomp**th.gm1) - 1.0 + self.dp2n
        # mpv.p2_nodes -= mpv.p2_nodes.mean()

        # pn = self.pn - self.pn.mean()
        # ps = self.ps - self.ps.mean()
        mpv.p2_nodes[...] = self.dp2n

    def walk_half_step(self):
        None

        # euler_forward_non_advective(Sol, mpv, elem, node, 0.5*dt, ud, th)

        # advect(Sol, flux, dt, elem, step%2, ud, th, mpv)

        # euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        # euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0, writer=None, label=str(label)+'_after_full_step')
