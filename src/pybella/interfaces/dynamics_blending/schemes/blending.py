import numpy as np
from scipy import signal


class Blend(object):
    """
    Class that takes care of the blending interface.
    """

    def __init__(self, ud):
        self.bb = False
        self.cb = ud.continuous_blending
        self.psinc_init = ud.no_of_pi_initial
        # self.psinc_trans = ud.no_of_pi_transition
        self.hydro_init = ud.no_of_hy_initial

        if self.psinc_init > 0 and self.cb:
            self.bb = True

        self.c_init = self.criterion_init
        # self.c_trans = self.criterion_trans

        self.fac = ud.Msq

    def criterion_init(self, step):
        return step == (self.psinc_init) and self.cb and self.bb

    # def criterion_trans(self, step):
    #     return step <= self.psinc_trans and self.cb and not self.bb

    def convert_p2n(self, p2n):
        ndim = p2n.ndim
        dp2n = p2n - p2n.mean()

        self.kernel = np.ones([2] * ndim)
        dp2c = signal.fftconvolve(dp2n, self.kernel, mode="valid") / self.kernel.sum()

        # self.dp2n = dp2n - dp2n.mean()
        # self.dp2c = dp2c - dp2c.mean()

        self.dp2n = dp2n
        self.dp2c = dp2c

        return dp2c

    def update_sol(self, mem, ud, sgn, label=None, writer=None):
        if writer != None:
            writer.populate(str(label) + "_before_blending", "dp2n", self.dp2n)
        if writer != None:
            writer.write_all(mem, str(label) + "_before_blending")

        sol = mem.sol
        npf = mem.npf
        th = mem.th

        if sgn == "bef":
            sign = -1.0
        elif sgn == "aft":
            sign = +1.0
        else:
            assert 0, "sgn == bef or sgn == aft"

        rho = np.copy(sol.rho)
        rhoY = np.copy(sol.rhoY)

        Y = rhoY / rho

        if ud.blending_mean == "rhoY":
            rhoYc = (rhoY**th.gm1 + sign * self.fac * self.dp2c) ** (th.gm1inv)
        elif ud.blending_mean == "1.0":
            rhoYc = (1.0 + sign * self.fac * self.dp2c) ** (th.gm1inv)

        alpha = rhoYc / sol.rhoY

        if ud.blending_conv == "rho":
            ### keep theta, convert rho
            sol.rho[...] = rho * alpha
            sol.rhoY[...] = sol.rho * Y

            rho_fac = sol.rho / rho
            sol.rhou[...] *= rho_fac
            sol.rhov[...] *= rho_fac
            sol.rhow[...] *= rho_fac
            sol.rhoX[...] *= rho_fac

        elif ud.blending_conv == "theta":
            ### keep rho, convert theta
            Yc = Y * alpha
            sol.rhoY[...] = rho * Yc
            sol.rhoX[...] = rho * (1.0 / Yc - npf.HydroState.S0.reshape(1, -1))
        else:
            assert 0, "ud.blending_conv undefined."

        if writer != None:
            writer.write_all(mem, str(label) + "_after_blending")

    def update_p2n(self, npf):
        npf.p2_nodes = self.dp2n
