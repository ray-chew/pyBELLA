import numpy as np

class NO_OF_RK_STAGES(object):
    NO_OF_RK_STAGES = 3

class TimeIntegratorParams(NO_OF_RK_STAGES):
    def __init__(self):
        self.dt_frac = 0
        self.flux_frac = np.zeros((self.NO_OF_RK_STAGES,2))
        self.update_frac = np.zeros((self.NO_OF_RK_STAGES))
        self.multiD_updt = False