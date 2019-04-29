import numpy as np
from input.enum_bdry import BdryType
from numerics_fundamentals.math_own import *

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
        if iny > 1:
            self.dy = (y1 - y0) / (iny - 1.)
        else:
            self.dy = 0.0
        if inz > 1:
            self.dz = (z1 - z0) / (inz - 1.)
        else:
            self.dz = 0.0
        
        assert self.dx > 0.0
        assert self.dy >= 0.0
        assert self.dz >= 0.0

        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.z0 = z0
        self.z1 = z1

        self.x = np.zeros((inx))
        self.y = np.zeros((iny))
        self.z = np.zeros((inz))

        for i in range(inx):
            self.x[i] = x0 + self.dx * i

        for i in range(iny):
            self.y[i] = y0 + self.dy * i

        for i in range(inz):
            self.z[i] = z0 + self.dz * i

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
            self.frotn = BdryType.TUNIX

big = 1.0

class ElemSpaceDiscr(object):
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
        if g.iny > 1:
            self.igy = self.ig[1] = 2
        else:
            self.igy = self.ig[1] = 0
        if g.inz > 1:
            self.igz = self.ig[2] = 2
        else:
            self.igz = self.ig[2] = 0

        self.icx = self.ic[0] = g.inx - 1 + 2 * self.igx
        if g.iny > 1:
            self.icy = self.ic[1] = g.iny - 1 + 2 * self.igy
        else:
            self.icy = 1
        if g.inz > 1:
            self.icz = self.ic[2] = g.inz - 1 + 2 * self.igz
        else:
            self.icz = 1

        self.nc = self.icx * self.icy * self.icz

        self.stride[0] = 1
        self.stride[1] = self.icx if g.iny > 1 else 0
        self.stride[2] = self.icx * self.icy if g.inz > 1 else 0

        self.ifx = self.icx + 1
        self.ify = self.icy + 1 if g.iny > 1 else 0
        self.ifz = self.icz + 1 if g.inz > 1 else 0

        self.nfx = self.ifx * self.icy * self.icz
        self.nfy = self.icx * self.ify * self.icz
        self.nfz = self.icx * self.icy * self.ifz

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

        self.x = np.zeros((self.icx))
        self.y = np.zeros((self.icy))
        self.z = np.zeros((self.icz))

        x0 = g.x0 - self.igx * self.dx + 0.5 * self.dx
        y0 = g.y0 - self.igy * self.dy + 0.5 * self.dy if self.icy > 1 else g.y0
        z0 = g.z0 - self.igz * self.dz + 0.5 * self.dz if self.icz > 1 else g.z0

        for i in range(self.icx):
            self.x[i] = x0 + self.dx * i
        for j in range(self.icy):
            self.y[i] = y0 + self.dy * j
        for k in range(self.icz):
            self.z[i] = z0 + self.dz * k

        self.left = g.left
        self.right = g.right
        self.bottom = g.bottom
        self.top = g.top
        self.back = g.back
        self.front = g.front

        self.scale_factor = 1.0
        