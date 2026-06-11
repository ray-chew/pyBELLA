"""JAX twin of :mod:`pybella.utils.operators.finite_difference`."""

import jax.numpy as jnp


def do_1d(field, spacing, axis=0):
    """
    Compute 1D finite difference along specified axis.

    Parameters
    ----------
    field : array
        Input field
    spacing : float
        Grid spacing
    axis : int, default=0
        Axis along which to compute difference

    Returns
    -------
    jax.Array
        Finite difference result
    """
    if field.ndim == 2:
        if axis == 0:
            return (field[1:, :] - field[:-1, :]) / spacing
        elif axis == 1:
            return (field[:, 1:] - field[:, :-1]) / spacing
        else:
            raise ValueError("axis must be 0 or 1 for 2D arrays")
    elif field.ndim == 3:
        if axis == 0:
            return (field[1:, :, :] - field[:-1, :, :]) / spacing
        elif axis == 1:
            return (field[:, 1:, :] - field[:, :-1, :]) / spacing
        elif axis == 2:
            return (field[:, :, 1:] - field[:, :, :-1]) / spacing
        else:
            raise ValueError("axis must be 0, 1, or 2 for 3D arrays")
    else:
        raise ValueError("field must be 2D or 3D array")
