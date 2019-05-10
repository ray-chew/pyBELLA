import numpy as np
from management.variable import States

class MPV(object):
    def __init__(self,elem,node,ud):
        sc = elem.sc
        sn = node.sc
        
        self.p0 = 1.0
        self.p00 = 1.0

        self.p2_cells = np.zeros((sc))
        self.dp2_cells = np.zeros((sc))
        self.p2_nodes = np.zeros((sn))
        self.p2_nodes0 = np.zeros((sn))
        self.dp2_nodes = np.zeros((sn))

        self.rhs = np.zeros((sc))
        self.diaginv = np.zeros((sc))
        self.wcenter = np.zeros((sn))
        self.wgrav = np.zeros((sn))
        self.wplus = np.zeros(([elem.ndim]+list(sn)))

        self.HydroState = States(sc,ud)
        self.HydroState_n = States(sn,ud)

        self.squeezer()

    def squeezer(self):
        for key, value in vars(self).items():
            if type(value) == np.ndarray:
                setattr(self,key,value.squeeze())