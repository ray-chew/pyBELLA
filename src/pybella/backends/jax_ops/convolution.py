"""JAX twin of :mod:`pybella.utils.operators.convolution`.

The kernel factories (``get_flux_kernels``, ``get_averaging_kernel``) and
``DIRECTION_CONFIG`` are static host-side data and are re-exported from the
numpy module unchanged.

Like the numba kernels, ``apply_convolution_kernel`` computes a "valid"
cross-correlation (no kernel flip). It is implemented as a shift-and-add
unrolled over the (tiny, static) kernel support in the same row-major order
as the numba loops, so results are bit-identical in x64.
"""

import numpy as np
import jax.numpy as jnp

from pybella.utils.operators.convolution import (
    DIRECTION_CONFIG,
    get_averaging_kernel,
    get_flux_kernels,
)

__all__ = [
    "DIRECTION_CONFIG",
    "get_averaging_kernel",
    "get_flux_kernels",
    "apply_convolution_kernel",
    "apply_directional_convolution",
]


def _correlate_valid(data, kernel):
    """Valid-mode cross-correlation, unrolled over the static kernel."""
    out_shape = tuple(d - k + 1 for d, k in zip(data.shape, kernel.shape))
    result = jnp.zeros(out_shape, dtype=jnp.float64)
    for kidx in np.ndindex(kernel.shape):
        sl = tuple(slice(i, i + n) for i, n in zip(kidx, out_shape))
        result = result + data[sl] * kernel[kidx]
    return result


def apply_convolution_kernel(
    data, kernel, normalize=True, axis_swap=None, use_numba=True
):
    """Apply convolution kernel with optional normalization and axis swapping.

    ``use_numba`` is accepted for signature parity with the numpy twin and
    ignored. ``kernel`` must be a concrete numpy array (it parameterizes the
    unrolled computation and is never traced).
    """
    kernel = np.asarray(kernel)
    if data.ndim not in (2, 3) or kernel.ndim != data.ndim:
        raise ValueError("data and kernel must both be 2D or 3D")

    result = _correlate_valid(jnp.asarray(data), kernel)

    if normalize:
        result = result / kernel.sum()

    if axis_swap is not None:
        result = jnp.moveaxis(result, axis_swap[0], axis_swap[1])

    return result


def apply_directional_convolution(
    data, kernel, direction, ndim, normalize=True, use_numba=True
):
    """Apply convolution kernel for a specific direction with axis swapping."""
    if direction not in DIRECTION_CONFIG[ndim]:
        raise ValueError(f"Direction '{direction}' not supported for {ndim}D")

    config = DIRECTION_CONFIG[ndim][direction]
    return apply_convolution_kernel(
        data,
        kernel,
        normalize=normalize,
        axis_swap=config["axis_swap"],
        use_numba=use_numba,
    )
