import numpy as np
import logging

from . import fields


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
        self.u = np.zeros((size), dtype=np.float64)
        self.v = np.zeros((size), dtype=np.float64)
        self.w = np.zeros((size), dtype=np.float64)
        self.X = np.zeros((size), dtype=np.float64)
        self.Y = np.zeros((size), dtype=np.float64)

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
        self._coriolis_cache = {}

    def get_recovery_objects(self, shape, ud):
        """Get cached recovery objects or create new ones."""

        cache_key = (tuple(shape), id(ud))
        if cache_key not in self._recovery_cache:
            logging.info("Cache: Creating new recovery objects with shape %s", shape)
            self._recovery_cache[cache_key] = {
                "Diffs": Characters(shape),
                "Ampls": Characters(shape),
                "Lefts": fields.States(shape, ud),
                "Rights": fields.States(shape, ud),
                "Slopes": Characters(shape),
            }

        cache_obj = self._recovery_cache[cache_key]

        return cache_obj

    def get_velocity_arrays(self, shape, dtype=np.float64):
        """Get cached velocity arrays (U, V, W) or create new ones."""

        cache_key = (tuple(shape), dtype)

        if cache_key not in self._velocity_cache:
            logging.info("Cache: Creating new velocity arrays with shape %s", shape)
            self._velocity_cache[cache_key] = {
                "U": np.zeros(shape, dtype=dtype),
                "V": np.zeros(shape, dtype=dtype),
                "W": np.zeros(shape, dtype=dtype),
            }

        cache_obj = self._velocity_cache[cache_key]

        return cache_obj

    def get_velocity_array_views(self, shape, dtype=np.float64):
        """Get views of cached velocity arrays for in-place operations."""
        cache_obj = self.get_velocity_arrays(shape, dtype)

        return cache_obj["U"], cache_obj["V"], cache_obj["W"]

    def get_coriolis_arrays(self, shape, dtype=np.float64):
        """Get cached Coriolis arrays (h11, h12, h13, h21, h22, h23, h31, h32, h33) or create new ones."""

        cache_key = (tuple(shape), dtype)

        if cache_key not in self._coriolis_cache:
            logging.info("Cache: Creating new Coriolis arrays with shape %s", shape)
            self._coriolis_cache[cache_key] = {
                "h11": np.zeros(shape, dtype=dtype),
                "h12": np.zeros(shape, dtype=dtype),
                "h13": np.zeros(shape, dtype=dtype),
                "h21": np.zeros(shape, dtype=dtype),
                "h22": np.zeros(shape, dtype=dtype),
                "h23": np.zeros(shape, dtype=dtype),
                "h31": np.zeros(shape, dtype=dtype),
                "h32": np.zeros(shape, dtype=dtype),
                "h33": np.zeros(shape, dtype=dtype),
                "denom": np.zeros(shape, dtype=dtype),
            }

        cache_obj = self._coriolis_cache[cache_key]

        return cache_obj

    def get_coriolis_array_views(self, shape, dtype=np.float64):
        """Get views of cached Coriolis arrays for in-place operations."""
        cache_obj = self.get_coriolis_arrays(shape, dtype)

        return (
            cache_obj["h11"],
            cache_obj["h12"],
            cache_obj["h13"],
            cache_obj["h21"],
            cache_obj["h22"],
            cache_obj["h23"],
            cache_obj["h31"],
            cache_obj["h32"],
            cache_obj["h33"],
            cache_obj["denom"],
        )

    def clear_all(self):
        """Clear all caches to free memory."""
        self._recovery_cache.clear()
        self._velocity_cache.clear()
