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

        self.u = np.zeros((sc))
        self.v = np.zeros((sc))
        self.w = np.zeros((sc))

        self.rhs = np.zeros((node.isc))
        self.wcenter = np.zeros((node.isc))
        self.wplus = np.zeros(([elem.ndim]+list(sc)))

        self.HydroState = States([sc[1]],ud)
        self.HydroState_n = States([sn[1]],ud)

        self.squeezer()

    def squeezer(self):
        for key, value in vars(self).items():
            if type(value) == np.ndarray:
                setattr(self,key,value.squeeze())


# def acoustic_order(ud,t,step):
#     if ud.is_compressible == 0:
#         return 2.0
#     elif ud.is_compressible == 1:
#         return 2.0
#     elif ud.is_compressible == -1:
#         current_transition_step = step - ud.no_of_pi_initial
#         # return np.linspace(,2.0,ud.no_of_pi_transition+2)[1:-1][current_transition_step]
#         return 2.0
#     else:
#         assert 0