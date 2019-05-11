import numpy as np
from inputs.enum_bdry import BdryType
from numerics_fundamentals.math_own import MIN_own

class Grid(object):
    def __init__(self, inx,iny,inz,x0,x1,y0,y1,z0,z1,left,right,bottom,top,back,front):
        assert inx > 1
        assert iny >= 1
        assert inz >= 1

        if inz > 1:
            self.ndim = 3
        elif iny > 1:
            self.ndim = 2
        else:
            self.ndim = 1
        
        self.inx = inx
        self.iny = iny
        self.inz = inz

        self.dx = (x1 - x0) / (inx - 1.)
        self.dy = (y1 - y0) / (iny - 1.) if iny > 1 else 0.0
        self.dz = (z1 - z0) / (inz - 1.) if inz > 1 else 0.0
        
        assert self.dx > 0.0
        assert self.dy >= 0.0
        assert self.dz >= 0.0

        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.z0 = z0
        self.z1 = z1

        self.x = x0 + self.dx * np.arange(inx)
        self.y = y0 + self.dy * np.arange(iny)
        self.z = z0 + self.dz * np.arange(inz)

        self.left = left
        self.right = right
        self.bottom = bottom
        self.top = top
        self.back = back
        self.front = front

        if iny == 1:
            self.bottom = BdryType.TUNIX
            self.top = BdryType.TUNIX

        if inz == 1:
            self.back = BdryType.TUNIX
            self.front = BdryType.TUNIX

big = 1.0

class SpaceDiscr(object):
    ig = np.zeros((3))
    ic = np.zeros((3))
    stride = np.zeros((3))
    dxyz = np.zeros((3))

    def __init__(self,g):
        assert g.inx > 1
        assert g.iny >= 1
        assert g.inz >= 1

        self.ndim = g.ndim
        self.normal = big

        self.igx = self.ig[0] = 2
        self.igy = self.ig[1] = 2 if g.iny > 1 else 0
        self.igz = self.ig[2] = 2 if g.inz > 1 else 0

        self.icx = self.ic[0] = g.inx - 1 + 2 * self.igx
        self.icy = self.ic[1] = g.iny - 1 + 2 * self.igy if g.iny > 1 else 1
        self.icz = self.ic[2] = g.inz - 1 + 2 * self.igz if g.inz > 1 else 1

        self.nc = self.icx * self.icy * self.icz
        # self.sc = (g.inx, g.iny, g.inz)
        self.sc = (self.icx, self.icy, self.icz)
        self.igs = [self.igx,self.igy,self.igz]

        self.stride[0] = 1
        self.stride[1] = self.icx if g.iny > 1 else 0
        self.stride[2] = self.icx * self.icy if g.inz > 1 else 0

        self.ifx = self.icx + 1
        self.ify = self.icy + 1 if g.iny > 1 else 0
        self.ifz = self.icz + 1 if g.inz > 1 else 0

        self.nfx = self.ifx * self.icy * self.icz
        self.nfy = self.icx * self.ify * self.icz
        self.nfz = self.icx * self.icy * self.ifz

        self.sfx = (self.ifx , self.icy , self.icz)
        self.sfy = (self.icx , self.ify , self.icz)
        self.sfz = (self.icx , self.icy , self.ifz)

        self.nf = self.nfx + self.nfy + self.nfz

        self.dx = g.dx
        self.dy = g.dy if self.icy > 1 else big
        self.dz = g.dz if self.icz > 1 else big
        
        self.dxyz[0] = self.dx
        self.dxyz[1] = self.dy
        self.dxyz[2] = self.dz

        assert self.dx > 0.0
        assert self.dy > 0.0
        assert self.dz > 0.0

        self.dxmin = MIN_own(self.dx, self.dy)
        self.dxmin = MIN_own(self.dxmin, self.dz)

        self.x = np.zeros((self.icx)).reshape(-1,1,1)
        self.y = np.zeros((self.icy)).reshape(1,-1,1)
        self.z = np.zeros((self.icz)).reshape(1,1,-1)

        self.left = g.left
        self.right = g.right
        self.bottom = g.bottom
        self.top = g.top
        self.back = g.back
        self.front = g.front

        self.scale_factor = 1.0
        
        self.inner_domain = np.empty((self.ndim),dtype=object)
        for dim in range(self.ndim):
            self.inner_domain[dim] = slice(self.igs[dim],-self.igs[dim])
        self.inner_domain = tuple(self.inner_domain)

        # save the indices of the ghost cells
        self.idx_ghost_min = np.empty((3), dtype=object)
        self.idx_ghost_max = np.empty((3), dtype=object)
        self.idx_ghost_min[0] = self.idx('i',0,self.igx)
        self.idx_ghost_max[0] = self.idx('i',-self.igx,)
        self.idx_ghost_min[1] = self.idx('j',0,self.igy)
        self.idx_ghost_max[1] = self.idx('j',-self.igy,)
        self.idx_ghost_min[2] = self.idx('k',0,self.igz)
        self.idx_ghost_max[2] = self.idx('k',-self.igz,)

        # save the indices of the inner boundary corresponding
        # to the number of ghost cells
        self.idx_inner_min = np.empty((3), dtype=object)
        self.idx_inner_max = np.empty((3), dtype=object)
        self.idx_inner_min[0] = self.idx('i',self.igx,self.igx+self.igx)
        self.idx_inner_max[0] = self.idx('i',-self.igx-self.igx,-self.igx)
        self.idx_inner_min[1] = self.idx('j',self.igy,self.igy+self.igy)
        self.idx_inner_max[1] = self.idx('j',-self.igy-self.igy,-self.igy)
        self.idx_inner_min[2] = self.idx('k',self.igz,self.igz+self.igz)
        self.idx_inner_max[2] = self.idx('k',-self.igz-self.igz,-self.igz)

        self.idx_inner_domain = np.empty((4), dtype=object)
        self.idx_inner_domain[0] = self.idx('i',self.igx,-self.igx)
        self.idx_inner_domain[1] = self.idx('j',self.igy,-self.igy)
        self.idx_inner_domain[2] = self.idx('k',self.igz,-self.igz)
        self.idx_inner_domain[3] = (self.idx_inner_domain[0][1], self.idx_inner_domain[1][0], self.idx_inner_domain[2][1])

    def idx(self,axs,start,stop=None):
        ijk = {
            'i' : -1,
            'j' : -2,
            'k' : -3,
        }
        # + ndim to account for python indexing; last axis = first index.
        axis = ijk[axs] + self.ndim

        indices = np.empty((self.ndim),dtype=object)
        indices[:] = slice(None)
        indices[axis] = slice(start,stop)
        return tuple(indices)
        

class ElemSpaceDiscr(SpaceDiscr):
    def __init__(self,g):
        super().__init__(g)
        x0 = g.x0 - self.igx * self.dx + 0.5 * self.dx
        y0 = g.y0 - self.igy * self.dy + 0.5 * self.dy if self.icy > 1 else g.y0
        z0 = g.z0 - self.igz * self.dz + 0.5 * self.dz if self.icz > 1 else g.z0

        self.x = x0 + self.dx * np.arange(self.icx)
        self.y = y0 + self.dy * np.arange(self.icy)
        self.z = z0 + self.dz * np.arange(self.icz)

        # for i in range(self.icx):
        #     self.x[i] = x0 + self.dx * i
        # for j in range(self.icy):
        #     self.y[i] = y0 + self.dy * j
        # for k in range(self.icz):
        #     self.z[i] = z0 + self.dz * k
        
class NodeSpaceDiscr(SpaceDiscr):
    def __init__(self,g):
        super().__init__(g)
        x0 = g.x0 - self.igx * self.dx
        y0 = g.y0 - self.igy * self.dy if self.icy > 1 else g.y0
        z0 = g.z0 - self.igz * self.dz if self.icz > 1 else g.z0

        self.icx += 1
        self.icy += 1 if g.iny > 1 else 0
        self.icz += 1 if g.inz > 1 else 0

        self.x = x0 + self.dx * np.arange(self.icx)
        self.y = y0 + self.dy * np.arange(self.icy)
        self.z = z0 + self.dz * np.arange(self.icz)

        self.sc = (self.icx , self.icy , self.icz)

        # for i in range(self.icx):
        #     self.x[i] = x0 + self.dx * i
        # for j in range(self.icy):
        #     self.y[i] = y0 + self.dy * j
        # for k in range(self.icz):
        #     self.z[i] = z0 + self.dz * k
