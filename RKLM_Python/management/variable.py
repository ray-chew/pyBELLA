import numpy as np
from scipy import signal

# equivalent to States_new
class Vars(object):
    """
    The data container for the solution state variables, i.e. `Sol`.

    """

    def __init__(self,size,ud):
        """
        Parameters
        ----------
        size : tuple
            Tuple containing the number of cells in the respective directions including ghost cells, e.g. `(48,48,10)` has 48 cells in the x and y-directions, and 10 cells in the z-directions
        ud : :class:`inputs.user_data.UserDataInit`
            Data container for the initial conditions

        Attributes
        ----------
        rho : ndarray(size)
        rhou : ndarray(size)
        rhov : ndarray(size)
        rhow : ndarray(size)
        rhoY : ndarray(size)
        rhoX : ndarray(size)

        Notes
        -----
        2. `rhoX` has to be extended by `ud.nspec` for moist process.
        
        """
        self.rho = np.zeros((size))
        self.rhou = np.zeros((size))
        self.rhov = np.zeros((size))
        self.rhow = np.zeros((size))
        self.rhoY = np.zeros((size))
        self.rhoX = np.zeros(([ud.nspec] + list(size)))
        self.squeezer()

    # will be a better way of doing this
    def squeezer(self):
        """
        Removes dimension of size 1. All arrays are initialised as 3D arrays, this function will remove the unnecessary dimensions.

        """
        for key, value in vars(self).items():
            setattr(self,key,value.squeeze())
    # method written for 2D

    def primitives(self,th):
        """
        Calculate the primitive quantities from the state variables and extend the data container to include these quantities.

        Parameters
        ----------
        th : :class:`physics.gas_dynamics.thermodynamic.ThermodynamicInit`
            Thermodynamic variables of the system

        Attributes
        ----------
        u : ndarray(size_of_rhou)
        v : ndarray(size_of_rhov)
        w : ndarray(size_of_rhow)
        Y : ndarray(size_of_rhoY)
        X : ndarray(size_of_rhoX)
        p : ndarray(size_of_rhoY)
     
        """
        nonzero_idx = np.nonzero(self.rho)

        self.u = np.zeros_like(self.rhou)
        self.v = np.zeros_like(self.rhov)
        self.w = np.zeros_like(self.rhow)
        self.Y = np.zeros_like(self.rhoY)
        self.X = np.zeros_like(self.rhoX)
        self.p = np.zeros_like(self.rhoY)

        # the non-zero indices are used here for the case where Lefts and
        # Rights in the HLLE Solver has one column that is zeros.
        self.u[nonzero_idx] = self.rhou[nonzero_idx] / self.rho[nonzero_idx]
        self.v[nonzero_idx] = self.rhov[nonzero_idx] / self.rho[nonzero_idx]
        self.w[nonzero_idx] = self.rhow[nonzero_idx] / self.rho[nonzero_idx]
        self.Y[nonzero_idx] = self.rhoY[nonzero_idx] / self.rho[nonzero_idx]
        self.X[nonzero_idx] = self.rhoX[nonzero_idx] / self.rho[nonzero_idx]
        self.p[nonzero_idx] = self.rhoY[nonzero_idx]**th.gamm

    def flip(self):
        """
        Flips the solution variables arrays for the advection routine. `rhou` and `rhov` are also flipped, i.e.::

            self.rhou, self.rhov = self.rhov, self.rhou    

        """
        for key, value in vars(self).items():
            setattr(self,key,value.T)

        self.rhou, self.rhov = self.rhov, self.rhou

    def flip_forward(self):
        for key, value in vars(self).items():
            setattr(self,key,np.moveaxis(value,0,-1))

    def flip_backward(self):
        for key, value in vars(self).items():
            setattr(self,key,np.moveaxis(value,-1,0))


    def mod_bg_wind(self, ud, fac):
        u0 = ud.u_wind_speed
        v0 = ud.v_wind_speed
        w0 = ud.w_wind_speed
        
        self.rhou = self.rhou + fac * u0 * self.rho
        self.rhov = self.rhov + fac * v0 * self.rho
        self.rhow = self.rhow + fac * w0 * self.rho
        
class States(Vars):
    """
    Data container for `Lefts` and `Rights` for the Riemann solver. Inherits the solution class :class:`management.variable.Vars`.

    """
    def __init__(self,size,ud):
        """
        Parameters
        ----------
        size : tuple
            Tuple containing the number of cells in the respective directions including ghost cells.
        ud : :class:`inputs.user_data.UserDataInit`
            Data container for the initial conditions

        Notes
        -----
        Many variables in this data container are no longer used and can be removed.

        """
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
        self.get_dSdy = self.get_dSdy
        self.get_S0c = self.get_S0c

    def get_dSdy(self, elem, node):
        if hasattr(self, 'dSdy'):
            return self.dSdy
        else:
            ndim = node.ndim
            dy = node.dy

            dSdy = self.S0
            dSdy = signal.convolve(dSdy,[1.,-1.],mode='valid') / dy

            for dim in range(0,ndim,2):
                dSdy = np.expand_dims(dSdy, dim)
                dSdy = np.repeat(dSdy, elem.sc[dim], axis=dim)

            self.dSdy = dSdy
            return dSdy

    def get_S0c(self, elem):
        if hasattr(self, 'S0c'):
            return self.S0c
        else:
            ndim = elem.ndim
            S0c = self.S0

            for dim in range(0,ndim,2):
                S0c = np.expand_dims(S0c, dim)
                S0c = np.repeat(S0c, elem.sc[dim], axis=dim)

            self.S0c = S0c
            return S0c


class Characters(object):
    """
    Data container for the slope and amplitude of the interpolation to the faces for the Riemann solver.

    """
    def __init__(self, size):
        """
        Parameters
        ----------
        size : tuple
            Tuple containing the number of cells in the respective directions including ghost cells.

        Attributes
        ----------
        u : ndarray(size)
        v : ndarray(size)
        w : ndarray(size)
        Y : ndarray(size)
        X : ndarray(size)
        plus : ndarray(size)
        minus : ndarray(size)
        entro : ndarray(size)

        """
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
        """
        Removes dimension of size 1. All arrays are initialised as 3D arrays, this function will remove the unnecessary dimensions.

        """
        for key, value in vars(self).items():
            setattr(self, key, value.squeeze())
