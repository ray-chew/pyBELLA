class ThermodynamicInit(object):
    """
    Data container for thermodynamical quantities.

    """
    def __init__(self, ud):
        """
        Parameters
        ----------
        ud : :class:`inputs.user_data.UserDataInit`
            Data container for the initial conditions
        """
        
        g = ud.gamm
        self.gamm = g
        self.gamminv = 1.0 / g
        self.gm1 = g - 1.0
        self.gm1inv = 1.0 / (g - 1.0)
        self.Gamma = (g - 1.0) / g
        self.Gammainv = g / (g - 1.0)

