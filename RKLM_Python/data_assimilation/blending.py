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

    def interface(self):
        None
        
    def convert_p2n(self, p2n):
        ndim =  p2n.ndim
        dp2n = p2n - p2n.mean()

        self.kernel = np.ones([2]*ndim)
        dp2c = signal.fftconvolve(dp2n,self.kernel, mode='valid') / self.kernel.sum()

        self.dp2n = dp2n - dp2n.mean()
        self.dp2c = dp2c - dp2c.mean()

        return dp2c

    def update_Sol(self, Sol, elem, node, th, ud, mpv, label=None,writer=None):
        rho = np.copy(Sol.rho)
        rhoY = np.copy(Sol.rhoY)

        Y = rhoY / rho

        rhoYc = (Sol.rhoY**th.gm1 + self.fac * self.dp2c)**(th.gm1inv)

        alpha = rhoYc / Sol.rhoY


        ### keep theta, convert rho
        # diff = alpha - Sol.rhoY
        # # locs = np.where(diff < 0.0005)

        # Sol.rho[...] = rho*alpha
        # Sol.rhoY[...] = (Sol.rho * Y)
        ### shift rhoY
        # Sol.rhoY[locs] -= diff[locs]

        # rho_fac = Sol.rho / rho
        # Sol.rhou[...] *= rho_fac
        # Sol.rhov[...] *= rho_fac
        # Sol.rhow[...] *= rho_fac
        # Sol.rhoX[...] *= rho_fac


        ### keep rho, convert theta
        # tol = 100.
        # diffY = alpha - Sol.rhoY
        # diffX = alpha - Sol.rhoX

        # locsY = np.where(np.abs(diffY) < tol)
        # locsX = np.where(np.abs(diffX) < tol)

        Yc = Y * alpha
        Sol.rhoY[...] = rho * Yc
        Sol.rhoX[...] = rho / Yc
        
        ### shift rhoY and rhoX
        # Sol.rhoY[locsY] -= diffY[locsY]
        # Sol.rhoX[locsX] -= diffX[locsX]


        writer.write_all(Sol,mpv,elem,node,th,str(label)+'_after_blending')
        # writer.populate(str(label)+'_before_blending', 'dp2n', self.dp2n)

    def update_p2n(self,Sol, mpv, node, th, ud):
        # rhoYncomp = np.zeros_like(self.dp2n)
        # rhoYncomp[1:-1,1:-1] = signal.fftconvolve(Sol.rhoY,self.kernel,mode='valid') / self.kernel.sum()
        # boundary.set_ghostnodes_p2(rhoYncomp,node,ud)
        # mpv.p2_nodes[...] = (rhoYncomp**th.gm1) - 1.0 + self.dp2n
        # mpv.p2_nodes -= mpv.p2_nodes.mean()

        mpv.p2_nodes = self.dp2n

        ### shift p2n
        # tol = 100.
        # diff = mpv.p2_nodes - 0.0
        # locs = np.where(np.abs(diff) < tol)
        # mpv.p2_nodes[locs] -= diff[locs]



    def walk_half_step(self, Sol, flux, mpv, elem, node, dt, ud):
        None
        # recompute_advective_fluxes(flux, Sol)

        # euler_forward_non_advective(Sol, mpv, elem, node, 0.5*dt, ud, th)

        # advect(Sol, flux, dt, elem, step%2, ud, th, mpv)

        # euler_backward_non_advective_expl_part(Sol, mpv, elem, 0.5*dt, ud, th)
        # euler_backward_non_advective_impl_part(Sol, mpv, elem, node, ud, th, t, 0.5*dt, 2.0, writer=None, label=str(label)+'_after_full_step')
