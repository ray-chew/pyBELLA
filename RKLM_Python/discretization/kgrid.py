import numpy as np
from inputs.enum_bdry import BdryType

class Grid(object):
    # def __init__(self, inx,iny,inz,x0,x1,y0,y1,z0,z1,left,right,bottom,top,back,front):
    """
    Base grid class, defines the extent of grid and the grid spacing.

    """
    def __init__(self, inx,iny,inz,x0,x1,y0,y1,z0,z1):
        """
        Parameters
        ----------
        inx : int
            Number of grid points in the x direction
        iny : int
            Number of grid points in the y direction
        inz : int
            Number of grid points in the z direction
        x0 : float
            Minimum extent in the x direction
        x1 : float
            Maximum extent in the x direction
        y0 : float
            Minimum extent in the y direction
        y1 : float
            Maximum extent in the y direction
        z0 : float
            Minimum extent in the z direction
        z1 : float
            Maximum extent in the z direction   
        
        """
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

        # self.left = left
        # self.right = right
        # self.bottom = bottom
        # self.top = top
        # self.back = back
        # self.front = front

big = 1.0
class SpaceDiscr(object):
    """
    For a given grid extent and number of grid-points, this class returns an equidistant discretised grid.

    """

    ig = np.zeros((3))
    ic = np.zeros((3))
    stride = np.zeros((3))
    dxyz = np.zeros((3))

    def __init__(self,g):
        """
        Parameters
        ----------
        g : :class:`discretization.kgrid.Grid`
            A grid class defining the grid extent and number of grid-points in each direction.

        """

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

        self.iicx = self.ic[0] = g.inx - 1
        self.iicy = self.ic[1] = g.iny - 1 if g.iny > 1 else 1
        self.iicz = self.ic[2] = g.inz - 1 if g.inz > 1 else 1

        self.nc = self.icx * self.icy * self.icz
        self.iisc = (self.iicx, self.iicy, self.iicz)
        self.isc = (self.iicx+self.igx, self.iicy+self.igy, self.iicz+self.igz)
        self.sc = (self.icx, self.icy, self.icz)
        self.igs = [self.igx,self.igy,self.igz]

        self.ifx = self.icx + 1
        self.ify = self.icy + 1 if g.iny > 1 else 0
        self.ifz = self.icz + 1 if g.inz > 1 else 0

        self.nfx = self.ifx * self.icy * self.icz
        self.nfy = self.icx * self.ify * self.icz
        self.nfz = self.icx * self.icy * self.ifz

        self.sfx = (self.icy , self.icz , self.ifx)
        self.sfy = (self.icz , self.icx , self.ify)
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
                   

        i1 = np.empty(self.ndim, dtype='object')
        i2 = np.empty(self.ndim, dtype='object')
        for dim in range(self.ndim):
            i1[dim] = slice(self.igs[dim]-1, (-self.igs[dim]+1))
            i2[dim] = slice(self.igs[dim],-self.igs[dim])
        self.i1 = tuple(i1)
        self.i2 = tuple(i2)

class ElemSpaceDiscr(SpaceDiscr):
    """
    Inherits the class :class:`discretization.kgrid.SpaceDiscr`. For a given grid extent and number of grid-points, this class returns an equidistant discretised cell-based grid.

    """

    def __init__(self,g, ud):
        """
        Parameters
        ----------
        g : :class:`discretization.kgrid.Grid`
            A grid class defining the grid extent and number of grid-points in each direction.

        """
        super().__init__(g)
        x0 = g.x0 - self.igx * self.dx + 0.5 * self.dx
        y0 = g.y0 - self.igy * self.dy + 0.5 * self.dy if self.icy > 1 else g.y0
        z0 = g.z0 - self.igz * self.dz + 0.5 * self.dz if self.icz > 1 else g.z0

        self.x = x0 + self.dx * np.arange(self.icx)
        self.y = y0 + self.dy * np.arange(self.icy)
        self.z = z0 + self.dz * np.arange(self.icz)

        self.get_p_indim(ud)

    def flip(self):
        self.dx, self.dy = self.dy, self.dx
        self.icx, self.icy = self.icy, self.icx
        self.ifx, self.ify = self.ify, self.ifx

    def get_p_indim(self, ud):
        ndim = self.ndim
        igs = self.igs
        p_isc = []
        pp1_isc = []

        eindim = np.empty((ndim),dtype='object')
        for dim in range(ndim):
            is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
            eindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)

            p_isc.append(self.isc[dim] + 2 * is_periodic)
            pp1_isc.append(self.isc[dim] + 2 * is_periodic + 2)

        self.periodic_indim = tuple(eindim)
        self.p_isc = tuple(p_isc)
        self.pp1_isc = tuple(pp1_isc)
        
        
class NodeSpaceDiscr(SpaceDiscr):
    """
    Inherits the class :class:`discretization.kgrid.SpaceDiscr`. For a given grid extent and number of grid-points, this class returns an equidistant discretised node-based grid.

    """
    def __init__(self,g, ud):
        """
        Parameters
        ----------
        g : :class:`discretization.kgrid.Grid`
            A grid class defining the grid extent and number of grid-points in each direction.

        """
        super().__init__(g)
        x0 = g.x0 - self.igx * self.dx
        y0 = g.y0 - self.igy * self.dy if self.icy > 1 else g.y0
        z0 = g.z0 - self.igz * self.dz if self.icz > 1 else g.z0

        self.icx += 1
        self.icy += 1 if g.iny > 1 else 0
        self.icz += 1 if g.inz > 1 else 0

        self.iicx += 1
        self.iicy += 1 if g.iny > 1 else 0
        self.iicz += 1 if g.inz > 1 else 0

        self.x = x0 + self.dx * np.arange(self.icx)
        self.y = y0 + self.dy * np.arange(self.icy)
        self.z = z0 + self.dz * np.arange(self.icz)

        self.iisc = (self.iicx , self.iicy , self.iicz)
        self.isc = (self.iicx+self.igx , self.iicy+self.igy , self.iicz+self.igz)
        self.sc = (self.icx , self.icy , self.icz)

        self.get_p_indim(ud)

    def flip(self):
        self.dx, self.dy = self.dy, self.dx
        self.icx, self.icy = self.icy, self.icx
        self.ifx, self.ify = self.ify, self.ifx

    def get_p_indim(self, ud):
        ndim = self.ndim
        igs = self.igs
        p_isc = []
        pp1_isc = []

        nindim = np.empty((ndim),dtype='object')
        for dim in range(ndim):
            is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
            nindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
            
            p_isc.append(self.isc[dim] + 2 * is_periodic)
            pp1_isc.append(self.isc[dim] + 2 * is_periodic + 2)

        self.periodic_indim = tuple(nindim)
        self.p_isc = tuple(p_isc)
        self.pp1_isc = tuple(pp1_isc)