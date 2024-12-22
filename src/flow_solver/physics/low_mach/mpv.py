import numpy as np

from ...utils import variable as var

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

        self.u = np.zeros((sc))
        self.v = np.zeros((sc))
        self.w = np.zeros((sc))

        self.rhs = np.zeros((node.isc))
        self.wcenter = np.zeros((node.isc))
        self.wplus = np.zeros(([elem.ndim]+list(sc)))

        self.HydroState = var.States([sc[1]],ud)
        self.HydroState_n = var.States([sn[1]],ud)

        self.squeezer()

    def squeezer(self):
        for key, value in vars(self).items():
            if type(value) == np.ndarray:
                setattr(self,key,value.squeeze())
