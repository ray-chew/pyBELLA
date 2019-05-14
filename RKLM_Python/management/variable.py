import numpy as np


# equivalent to States_new
class Vars(object):
    def __init__(self,size,ud):
        self.rho = np.zeros((size))
        self.rhou = np.zeros((size))
        self.rhov = np.zeros((size))
        self.rhow = np.zeros((size))
        self.rhoe = np.zeros((size))
        self.rhoY = np.zeros((size))
        self.rhoX = np.zeros(([ud.nspec] + list(size)))
        self.squeezer()

    # will be a better way of doing this
    def squeezer(self):
        for key, value in vars(self).items():
            setattr(self,key,value.squeeze())
    # method written for 2D

class States(Vars):
    def __init__(self,size,ud):
        super().__init__(size,ud)
        self.u = np.zeros((size))
        self.v = np.zeros((size))
        self.w = np.zeros((size))
        self.q = np.zeros((size))
        self.p = np.zeros((size))
        self.c = np.zeros((size))
        self.entro = np.zeros((size))
        self.H = np.zeros((size))
        self.Y = np.zeros((size))
        self.X = np.zeros(([ud.nspec] + list(size)))

        self.p0 = np.zeros((size))
        self.p20 = np.zeros((size))
        self.rho0 = np.zeros((size))
        self.S0 = np.zeros((size))
        self.S10 = np.zeros((size))
        self.pi0 = np.zeros((size))
        self.rhoY0 = np.zeros((size))
        self.Y0 = np.zeros((size))

        self.squeezer()

# class StatesSmall(States):
#     def __init__(self,size,ud):
#         super().__init__(size,ud)


#         self.squeezer()
        
