import numpy as np
from management.variable import States

class MPV(object):
    def __init__(self,elem,node,ud):
        self.nc = elem.nc
        self.nn = node.nc

        self.p0 = 1.0
        self.p00 = 1.0

        self.p2_cells = np.zeros((self.nc))
        self.dp2_cells = np.zeros((self.nc))
        self.p2_nodes = np.zeros((self.nn))
        self.p2_nodes0 = np.zeros((self.nn))
        self.dp2_nodes = np.zeros((self.nn))

        self.rhs = np.zeros((self.nc))
        self.diaginv = np.zeros((self.nc))
        self.wcenter = np.zeros((self.nn))
        self.wgrav = np.zeros((self.nn))
        self.wplus = np.zeros((elem.ndim,self.nn))

        self.HydroState = States(self.nc,ud)
        self.HydroState_n = States(self.nn,ud)

