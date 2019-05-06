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
        self.rhoX = np.zeros((ud.nspec, size))

    # method written for 2D
    def set_periodic(self,elem):
        igx = elem.igx

        # can be rewritten better in terms of a function.
        self.rho[:,:igx], self.rho[:,-igx:] = self.rho[:,-igx-igx:-igx], self.rho[:,igx+1:igx+igx+1].copy()
        self.rhou[:,:igx], self.rhou[:,-igx:] = self.rhou[:,-igx-igx:-igx], self.rhou[:,igx+1:igx+igx+1].copy()
        self.rhov[:,:igx], self.rhov[:,-igx:] = self.rhov[:,-igx-igx:-igx], self.rhov[:,igx+1:igx+igx+1].copy()
        self.rhow[:,:igx], self.rhow[:,-igx:] = self.rhow[:,-igx-igx:-igx], self.rhow[:,igx+1:igx+igx+1].copy()
        self.rhoe[:,:igx], self.rhoe[:,-igx:] = self.rhoe[:,-igx-igx:-igx], self.rhoe[:,igx+1:igx+igx+1].copy()
        self.rhoY[:,:igx], self.rhoY[:,-igx:] = self.rhoY[:,-igx-igx:-igx], self.rhoY[:,igx+1:igx+igx+1].copy()
        self.rhoX[:,:,0:igx], self.rhoX[:,:,-igx:] = self.rhoX[:,:,-igx-igx:-igx], self.rhoX[:,:,igx+1:igx+igx+1].copy()

    def set_wall(self,elem):
        igx = elem.igx

        # again, can be written in a swap_columns function and a for loop... but first this.
        self.rho[:,:igx], self.rho[:,-igx:] = self.rho[:,igx:igx+igx][:,::-1], self.rho[:,-igx-igx:-igx][:,::-1].copy()
        self.rhou[:,:igx], self.rhou[:,-igx:] = -1.*self.rhou[:,igx:igx+igx][:,::-1], -1.*self.rhou[:,-igx-igx:-igx][:,::-1].copy()
        self.rhov[:,:igx], self.rhov[:,-igx:] = self.rhov[:,igx:igx+igx][:,::-1], self.rhov[:,-igx-igx:-igx][:,::-1].copy()
        self.rhow[:,:igx], self.rhow[:,-igx:] = self.rhow[:,igx:igx+igx][:,::-1], self.rhow[:,-igx-igx:-igx][:,::-1].copy()
        self.rhoe[:,:igx], self.rhoe[:,-igx:] = self.rhoe[:,igx:igx+igx][:,::-1], self.rhoe[:,-igx-igx:-igx][:,::-1].copy()
        self.rhoY[:,:igx], self.rhoY[:,-igx:] = self.rhoY[:,igx:igx+igx][:,::-1], self.rhoY[:,-igx-igx:-igx][:,::-1].copy()


# equivalent to ConsVars_new
class StatesSmall(Vars):
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
        self.X = np.zeros((ud.nspec, size))

        self.p0 = np.zeros((size))
        self.p20 = np.zeros((size))
        self.rho0 = np.zeros((size))
        self.S0 = np.zeros((size))
        self.S10 = np.zeros((size))
        self.Y0 = np.zeros((size))

class States(StatesSmall):
    def __init__(self,size,ud):
        super().__init__(size,ud)
        self.pi0 = np.zeros((size))
        self.rhoY0 = np.zeros((size))

        
