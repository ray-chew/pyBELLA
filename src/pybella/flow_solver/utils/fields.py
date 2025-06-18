import numpy as np
import scipy as sp
import logging


class CellSolField(object):
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
        with np.errstate(divide="ignore", invalid="ignore"):
            # Direct division without nonzero indexing
            # We know that when this method is called in recovery, we always have one column of zeroes in self.rho.
            self.u[...] = self.rhou / self.rho
            self.v[...] = self.rhov / self.rho
            self.w[...] = self.rhow / self.rho
            self.Y[...] = self.rhoY / self.rho
            self.X[...] = self.rhoX / self.rho
            self.p[...] = self.rhoY**th.gamm

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


class States(CellSolField):
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
        if not self.init_dSdy:
            logging.info("Computing dSdy")
            self.dSdy = sp.signal.convolve(self.S0, [1.0, -1.0], mode="valid") / node.dy

            for dim in range(0, node.ndim, 2):
                self.dSdy = np.expand_dims(self.dSdy, dim)
                self.dSdy = np.repeat(self.dSdy, elem.sc[dim], axis=dim)

            self.init_dSdy = True

        return self.dSdy

    def get_S0c(self, elem):
        if not self.init_S0c:
            logging.info("Computing S0c")
            S0c_result = self.S0

            for dim in range(0, elem.ndim, 2):
                S0c_result = np.expand_dims(S0c_result, dim)
                S0c_result = np.repeat(S0c_result, elem.sc[dim], axis=dim)

            self.S0c = S0c_result
            self.init_S0c = True

        return self.S0c


class NodePressureField(object):
    def __init__(self, elem, node, ud):
        sc = elem.sc
        sn = node.sc

        self.p0 = 1.0
        self.p00 = 1.0

        self.p2_cells = np.zeros((sc))
        self.dp2_cells = np.zeros((sc))
        self.p2_nodes = np.zeros((sn))
        self.p2_nodes0 = np.zeros((sn))
        self.dp2_nodes = np.zeros((sn))

        self.u = np.zeros((sc))
        self.v = np.zeros((sc))
        self.w = np.zeros((sc))

        self.rhs = np.zeros((node.isc))
        self.wcenter = np.zeros((node.isc))
        self.wplus = np.zeros(([elem.ndim] + list(sc)))

        self.HydroState = States([sc[1]], ud)
        self.HydroState_n = States([sn[1]], ud)

        self.squeezer()

    def squeezer(self):
        for key, value in vars(self).items():
            if type(value) == np.ndarray:
                setattr(self, key, value.squeeze())
