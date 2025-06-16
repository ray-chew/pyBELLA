import numpy as np
import scipy as sp
import logging

class Vars(object):
    """
    The data container for the solution state variables, i.e. `Sol`.

    """

    def __init__(self, size, ud):
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

        self.u = np.zeros((size))
        self.v = np.zeros((size))
        self.w = np.zeros((size))
        self.Y = np.zeros((size))
        self.X = np.zeros(([ud.nspec] + list(size)))
        self.p = np.zeros((size))

        self.squeezer()

    # will be a better way of doing this
    def squeezer(self):
        """
        Removes dimension of size 1. All arrays are initialised as 3D arrays, this function will remove the unnecessary dimensions.

        """
        for key, value in vars(self).items():
            setattr(self, key, value.squeeze())

    def primitives(self, th):
        """
        Calculate the primitive quantities from the state variables and extend the data container to include these quantities.

        Parameters
        ----------
        th : :class:`physics.gas_dynamics.thermodynamic.init`
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
        with np.errstate(divide='ignore', invalid='ignore'):
            # Direct division without nonzero indexing
            # We know that when this method is called in recovery, we always have one column of zeroes in self.rho.
            self.u[...] = self.rhou / self.rho
            self.v[...] = self.rhov / self.rho
            self.w[...] = self.rhow / self.rho
            self.Y[...] = self.rhoY / self.rho
            self.X[...] = self.rhoX / self.rho
            self.p[...] = self.rhoY ** th.gamm

    def flip(self):
        """
        Flips the solution variables arrays for the advection routine. `rhou` and `rhov` are also flipped, i.e.::

            self.rhou, self.rhov = self.rhov, self.rhou

        """
        for key, value in vars(self).items():
            setattr(self, key, value.T)

        self.rhou, self.rhov = self.rhov, self.rhou

    def flip_forward(self):
        for key, value in vars(self).items():
            setattr(self, key, np.moveaxis(value, 0, -1))

    def flip_backward(self):
        for key, value in vars(self).items():
            setattr(self, key, np.moveaxis(value, -1, 0))

    def mod_bg_wind(self, ud, fac):
        u0 = ud.u_wind_speed
        v0 = ud.v_wind_speed
        w0 = ud.w_wind_speed

        self.rhou[...] = self.rhou + fac * u0 * self.rho
        self.rhov[...] = self.rhov + fac * v0 * self.rho
        self.rhow[...] = self.rhow + fac * w0 * self.rho


class States(Vars):
    """
    Data container for `Lefts` and `Rights` for the Riemann solver. Inherits the solution class :class:`management.variable.Vars`.

    """

    def __init__(self, size, ud):
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
        super().__init__(size, ud)

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

        self.init_dSdy = False
        self.init_S0c = False

    def get_dSdy(self, elem, node):
        if self.init_dSdy:
            return self.dSdy
        else:
            ndim = node.ndim
            dy = node.dy

            dSdy = self.S0
            dSdy = sp.signal.convolve(dSdy, [1.0, -1.0], mode="valid") / dy

            for dim in range(0, ndim, 2):
                dSdy = np.expand_dims(dSdy, dim)
                dSdy = np.repeat(dSdy, elem.sc[dim], axis=dim)

            self.dSdy = dSdy
            self.init_dSdy = True
            return dSdy

    def get_S0c(self, elem):
        if self.init_S0c:
            return self.S0c
        else:
            ndim = elem.ndim
            S0c = self.S0

            for dim in range(0, ndim, 2):
                S0c = np.expand_dims(S0c, dim)
                S0c = np.repeat(S0c, elem.sc[dim], axis=dim)

            self.S0c = S0c
            self.init_S0c = True
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

        self.squeezer()

    def squeezer(self):
        """
        Removes dimension of size 1. All arrays are initialised as 3D arrays, this function will remove the unnecessary dimensions.

        """
        for key, value in vars(self).items():
            setattr(self, key, value.squeeze())


class FlowSolverCache:
    """Cache for flow solver specific computations."""
    
    def __init__(self):
        self._recovery_cache = {}
        self._velocity_cache = {}
    
    def get_recovery_objects(self, shape, ud):
        """Get cached recovery objects or create new ones."""
        
        cache_key = (tuple(shape), id(ud))
        if cache_key not in self._recovery_cache:
            logging.info("Cache: Creating new recovery objects with shape %s", shape)
            self._recovery_cache[cache_key] = {
                'Diffs': Characters(shape),
                'Ampls': Characters(shape),
                'Lefts': States(shape, ud),
                'Rights': States(shape, ud),
                'Slopes': Characters(shape),   
            }
        
        # Reset objects if they have reset methods
        cache_obj = self._recovery_cache[cache_key]
        # for obj in cache_obj.values():
        #     if hasattr(obj, 'zero'):
        #         obj.zero()
        
        return cache_obj
    
    def get_velocity_arrays(self, shape, dtype=np.float64):
        """Get cached velocity arrays (U, V, W) or create new ones."""
        
        cache_key = (tuple(shape), dtype)
        
        if cache_key not in self._velocity_cache:
            logging.info("Cache: Creating new velocity arrays with shape %s", shape)
            self._velocity_cache[cache_key] = {
                'U': np.zeros(shape, dtype=dtype),
                'V': np.zeros(shape, dtype=dtype),
                'W': np.zeros(shape, dtype=dtype)
            }
        
        # Clear arrays for reuse
        cache_obj = self._velocity_cache[cache_key]
        # for arr in cache_obj.values():
        #     arr.fill(0.0)
        
        return cache_obj
    
    def get_velocity_array_views(self, shape, dtype=np.float64):
        """Get views of cached velocity arrays for in-place operations."""
        cache_obj = self.get_velocity_arrays(shape, dtype)
        
        return cache_obj['U'], cache_obj['V'], cache_obj['W']
    
    def clear_velocity_cache(self):
        """Clear velocity cache to free memory."""
        self._velocity_cache.clear()
    
    def clear_all(self):
        """Clear all caches to free memory."""
        self._recovery_cache.clear()
        self._velocity_cache.clear()

