import numpy as np
from input.enum_bdry import BdryType

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


