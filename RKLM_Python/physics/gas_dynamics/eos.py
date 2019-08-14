import numpy as np

def nonhydrostasy(ud,t):
    if ud.is_nonhydrostatic == 0:
        return 0.0
    elif ud.is_nonhydrostatic == 1:
        return 1.0
    elif ud.is_nonhydrostatic == -1:
        a = 12.5
        b = 1.0 / 24.0
        c = np.min(1.0, np.max(0.0, b * (t-a)))
        return c
    else:
        assert 0

def compressibility(ud,t):
    if ud.is_compressible == 0:
        return 0,0
    elif ud.is_compressible == 1:
        return 1.0
    elif ud.is_compressible == -1:
        dtloc = 12.0
        a = 0.0 * dtloc
        b = 1.0 / dtloc / 20.0
        tau = np.min(1.0, np.max(0.0 * b * (t-a)))
        c = 0.5 * (1.0 - np.cos(np.pi * tau))
        return c
    else:
        assert 0


def rhoe(rho,u,v,w,p,ud,th):
    Msq = ud.compressibility * ud.Msq
    gm1inv = th.gm1inv

    return p * gm1inv + 0.5 * Msq * rho * (u**2 + v**2 + w**2)
