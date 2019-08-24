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

    def primitives(self,th):
        nonzero_idx = np.nonzero(self.rho)

        self.u = np.zeros_like(self.rhou)
        self.v = np.zeros_like(self.rhov)
        self.w = np.zeros_like(self.rhow)
        self.Y = np.zeros_like(self.rhoY)
        self.X = np.zeros_like(self.rhoX)
        self.p = np.zeros_like(self.rhoY)

        self.u[nonzero_idx] = self.rhou[nonzero_idx] / self.rho[nonzero_idx]
        self.v[nonzero_idx] = self.rhov[nonzero_idx] / self.rho[nonzero_idx]
        self.w[nonzero_idx] = self.rhow[nonzero_idx] / self.rho[nonzero_idx]
        self.Y[nonzero_idx] = self.rhoY[nonzero_idx] / self.rho[nonzero_idx]
        self.X[nonzero_idx] = self.rhoX[nonzero_idx] / self.rho[nonzero_idx]
        self.p[nonzero_idx] = self.rhoY[nonzero_idx]**th.gamm

    def flip(self):
        for key, value in vars(self).items():
            setattr(self,key,value.T)

        self.rhou, self.rhov = self.rhov, self.rhou
        


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

class Characters(object):
    def __init__(self, size):
        self.u = np.zeros((size))
        self.v = np.zeros((size))
        self.w = np.zeros((size))
        self.X = np.zeros((size))
        self.Y = np.zeros((size))

        self.plus = np.zeros((size))
        self.minus = np.zeros((size))
        self.entro = np.zeros((size))

        self.squeezer()

    def squeezer(self):
        for key, value in vars(self).items():
            setattr(self, key, value.squeeze())

    def change_dir(self):
        for key, value in vars(self).items():
            setattr(self, key, -1. * value)
        
# class StatesSmall(States):
#     def __init__(self,size,ud):
#         super().__init__(size,ud)


#         self.squeezer()
        
