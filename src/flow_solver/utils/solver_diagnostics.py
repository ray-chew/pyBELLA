import numpy as np
import scipy as sp

def get_p_from_pressure_related_fields(mem, ud, psinc=False):
    p2n = mem.mpv.p2_nodes
    rhoY = mem.sol.rhoY
    dp2n = (p2n - p2n.mean()) * ud.Msq

    th = mem.th

    kernel = np.ones((2,2))
    dp2c = sp.signal.fftconvolve(dp2n, kernel, mode='valid') / kernel.sum()

    if psinc:
        P0 = (rhoY**(th.gamm-1.0) + dp2c)**(1.0/(th.gamm-1.0))
        p = P0**(th.gamm)
    else:
        P0 = (rhoY**(th.gamm-1.0) - dp2c)**(1.0/(th.gamm-1.0))
        p = rhoY**(th.gamm) - P0**(th.gamm)

    # non-dimensionalised p/p_ref
    return p