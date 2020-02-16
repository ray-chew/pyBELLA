class ThemodynamicInit(object):
    def __init__(self, ud):
        g = ud.gamm
        self.gamm = g
        self.gamminv = 1.0 / g
        self.gm1 = g - 1.0
        self.gm1inv = 1.0 / (g - 1.0)
        self.Gamma = (g - 1.0) / g
        self.Gammainv = g / (g - 1.0)

        self.Rg_over_Rv = ud.Rg_over_Rv
        self.Q = ud.Q

