import numpy as np
from enum import Enum
from physics.gas_dynamics.explicit import TimeIntegratorParams

class Test(object):
    def __init__(self,**kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

class TravellingVortex3D_48(object):
    h_ref = 10000
    R_gas = 287.4
    R_vap = 461.0

    def __init__(self):
        self.h_ref = self.h_ref
        self.R_gas = self.R_gas
        self.R_vap = self.R_vap
        self.Rg_over_Rv = self.R_gas / self.R_vap
    
class MolecularTransport(Enum):
    FULL_MOLECULAR_TRANSPORT = 0
    STRAKA_DIFFUSION_MODEL = 1
    NO_MOLECULAR_TRANSPORT = 2

variables = vars(TravellingVortex3D_48())
D = Test(**variables)
print(vars(D))

MTs = MolecularTransport.FULL_MOLECULAR_TRANSPORT
print(MTs == MolecularTransport.NO_MOLECULAR_TRANSPORT)

NO_OF_RK_STAGES = 3

class Ud(object):
    def __init__(self):
        self.tips = TimeIntegratorParams()
        print(self.tips.flux_frac)
        self.tips.flux_frac = 100
        print(self.tips.flux_frac)

ud = Ud()
def change_ud(ud):
    ud.tips.flux_frac = 5
change_ud(ud)
print(ud.tips.flux_frac)