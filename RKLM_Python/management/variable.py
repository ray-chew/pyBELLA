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

    def set_boundary(self,pads,btype,idx):
        self.rho[...] = np.pad(self.rho[idx],pads,btype)
        if btype == 'symmetric':
            self.rhou[...] = np.pad(self.rhou[idx],pads,negative_symmetric)
        else:
            self.rhou[...] = np.pad(self.rhou[idx],pads,btype)
        self.rhov[...] = np.pad(self.rhov[idx],pads,btype)
        self.rhow[...] = np.pad(self.rhow[idx],pads,btype)
        self.rhoe[...] = np.pad(self.rhoe[idx],pads,btype)
        self.rhoY[...] = np.pad(self.rhoY[idx],pads,btype)
        self.rhoX[...] = np.pad(self.rhoX[idx],pads,btype)

def negative_symmetric(vector,pad_width,iaxis,kwargs=None):
    if pad_width[1] > 0:
        sign = -1
        vector[:pad_width[0]] = sign * vector[pad_width[0]:2*pad_width[0]][::-1]
        vector[-pad_width[1]:] = sign * vector[-2*pad_width[1]:-pad_width[1]][::-1]
        return vector
    else:
        return vector

    #     # slice_along = {
    #     #     '0' : 1,
    #     #     '1' : 0,
    #     #     '2' : 2,
    #     # }
    #     # axs = slice_along[axs]


    #     # how to loop over attributes of a class?
    #     # possible to flip the arrays and use ellipses

    #     min_ghost = elem.idx_ghost_min[axs]
    #     max_ghost = elem.idx_ghost_max[axs]
    #     min_inner = elem.idx_inner_min[axs]
    #     max_inner = elem.idx_inner_max[axs]

    #     self.rho[min_ghost], self.rho[max_ghost] = self.rho[max_inner], self.rho[min_inner].copy()
    #     self.rhou[min_ghost], self.rhou[max_ghost] = self.rhou[max_inner], self.rhou[min_inner].copy()
    #     self.rhov[min_ghost], self.rhov[max_ghost] = self.rhov[max_inner], self.rhov[min_inner].copy()
    #     self.rhow[min_ghost], self.rhow[max_ghost] = self.rhow[max_inner], self.rhow[min_inner].copy()
    #     self.rhoe[min_ghost], self.rhoe[max_ghost] = self.rhoe[max_inner], self.rhoe[min_inner].copy()
    #     self.rhoY[min_ghost], self.rhoY[max_ghost] = self.rhoY[max_inner], self.rhoY[min_inner].copy()
    #     self.rhoX[:,:igx], self.rhoX[:,:,-igx:] = self.rhoX[:,max_inner], self.rhoX[:,min_inner].copy()

    # def set_wall(self,elem):
    #     igx = elem.igx

    #     min_ghost = elem.idx(axs,0,igx)
    #     max_ghost = elem.idx(axs,-igx,)
    #     min_inner = elem.idx(axs,igx,igx+igx+1)
    #     max_inner = elem.idx(axs,-igx-igx,-igx)

    #     # again, can be written in a swap_columns function and a for loop... but first this.
    #     self.rho[min_ghost], self.rho[max_ghost] = self.rho[min_inner][:,::-1], self.rho[max_inner][:,::-1].copy()
    #     self.rhou[min_ghost], self.rhou[max_ghost] = -1.*self.rhou[min_inner][:,::-1], -1.*self.rhou[max_inner][:,::-1].copy()
    #     self.rhov[min_ghost], self.rhov[max_ghost] = self.rhov[min_inner][:,::-1], self.rhov[max_inner][:,::-1].copy()
    #     self.rhow[min_ghost], self.rhow[max_ghost] = self.rhow[min_inner][:,::-1], self.rhow[max_inner][:,::-1].copy()
    #     self.rhoe[min_ghost], self.rhoe[max_ghost] = self.rhoe[min_inner][:,::-1], self.rhoe[max_inner][:,::-1].copy()
    #     self.rhoY[min_ghost], self.rhoY[max_ghost] = self.rhoY[min_inner][:,::-1], self.rhoY[max_inner][:,::-1].copy()


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
        
