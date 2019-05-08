import numpy as np
from management.variable import States

class MPV(object):
    def __init__(self,elem,node,ud):
        self.nc = elem.nc
        self.nn = node.nc

        self.sc = elem.sc
        self.sn = node.sc

        self.p0 = 1.0
        self.p00 = 1.0

        self.p2_cells = np.zeros((self.sc))
        self.dp2_cells = np.zeros((self.sc))
        self.p2_nodes = np.zeros((self.sn))
        self.p2_nodes0 = np.zeros((self.sn))
        self.dp2_nodes = np.zeros((self.sn))

        self.rhs = np.zeros((self.sc))
        self.diaginv = np.zeros((self.sc))
        self.wcenter = np.zeros((self.sn))
        self.wgrav = np.zeros((self.sn))
        self.wplus = np.zeros(([elem.ndim]+list(self.sn)))

        self.HydroState = States(self.sc,ud)
        self.HydroState_n = States(self.sn,ud)