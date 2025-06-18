import numpy as np
import numba as nb
import scipy as sp
import functools

# Configuration for directional convolutions
DIRECTION_CONFIG = {
    2: {"u": {"axis_swap": (0, -1)}, "v": {"axis_swap": None}},
    3: {
        "u": {"axis_swap": (0, -1)},
        "v": {"axis_swap": (-1, 0)},
        "w": {"axis_swap": None},
    },
}

@functools.lru_cache(maxsize=2)
def get_flux_kernels(ndim):
    """Create convolution kernels for advective flux computation.

    Parameters
    ----------
    ndim : int
        Number of dimensions (2 or 3)

    Returns
    -------
    dict
        Dictionary containing kernels for each direction

    Notes
    -----
    Results are cached since kernels don't change during computation.
    """
    if ndim == 2:
        kernel_u = np.array([[0.5, 1.0, 0.5], [0.5, 1.0, 0.5]])
        return {"u": kernel_u, "v": kernel_u.T}
    elif ndim == 3:
        kernel_u = np.array(
            [[[1, 2, 1], [2, 4, 2], [1, 2, 1]], [[1, 2, 1], [2, 4, 2], [1, 2, 1]]]
        )
        return {
            "u": kernel_u,
            "v": np.swapaxes(kernel_u, 1, 0),
            "w": np.swapaxes(kernel_u, 2, 0),
        }
    else:
        raise ValueError(f"Unsupported dimension: {ndim}")


@functools.lru_cache(maxsize=4)
def get_averaging_kernel(ndim, width=3, normalize=True):
    """
    Create a generic averaging kernel for arbitrary dimensions.

    Parameters
    ----------
    ndim : int
        Number of dimensions (e.g., 2 or 3)
    width : int, default=3
        Size of the kernel along each axis (can be even or odd)
    normalize : bool, default=True
        Whether to normalize the kernel to sum to 1

    Returns
    -------
    np.ndarray
        Averaging kernel of shape (width,) * ndim

    Notes
    -----
    - Odd widths result in centered kernels.
    - Even widths are useful for staggered/grid-face averaging.
    """
    shape = (width,) * ndim
    kernel = np.ones(shape, dtype=np.float64)

    if normalize:
        kernel /= kernel.size

    return kernel


@nb.njit(cache=True)
def _convolve_2d(data, kernel):
    """Numba-compiled 2D convolution for better performance."""
    data_h, data_w = data.shape
    kernel_h, kernel_w = kernel.shape

    result_h = data_h - kernel_h + 1
    result_w = data_w - kernel_w + 1
    result = np.zeros((result_h, result_w))

    for i in range(result_h):
        for j in range(result_w):
            for ki in range(kernel_h):
                for kj in range(kernel_w):
                    result[i, j] += data[i + ki, j + kj] * kernel[ki, kj]

    return result


@nb.njit(cache=True)
def _convolve_3d(data, kernel):
    """Numba-compiled 3D convolution for better performance."""
    data_d, data_h, data_w = data.shape
    kernel_d, kernel_h, kernel_w = kernel.shape

    result_d = data_d - kernel_d + 1
    result_h = data_h - kernel_h + 1
    result_w = data_w - kernel_w + 1
    result = np.zeros((result_d, result_h, result_w))

    for i in range(result_d):
        for j in range(result_h):
            for k in range(result_w):
                for ki in range(kernel_d):
                    for kj in range(kernel_h):
                        for kk in range(kernel_w):
                            result[i, j, k] += (
                                data[i + ki, j + kj, k + kk] * kernel[ki, kj, kk]
                            )

    return result


def apply_convolution_kernel(
    data, kernel, normalize=True, axis_swap=None, use_numba=True
):
    """Apply convolution kernel with optional normalization and axis swapping.

    Parameters
    ----------
    data : np.ndarray
        Input data array
    kernel : np.ndarray
        Convolution kernel
    normalize : bool, default=True
        Whether to normalize by kernel sum
    axis_swap : tuple or None, default=None
        Tuple of (from_axis, to_axis) for np.moveaxis
    use_numba : bool, default=True
        Whether to use Numba-compiled convolution (faster for repeated calls)

    Returns
    -------
    np.ndarray
        Convolved result

    Notes
    -----
    For large arrays or single calls, scipy.signal.fftconvolve might be faster.
    For repeated calls on smaller arrays, Numba convolution is typically faster.
    """
    if use_numba and data.ndim in (2, 3):
        if data.ndim == 2:
            result = _convolve_2d(data, kernel)
        else:  # 3D
            result = _convolve_3d(data, kernel)
    else:
        # Fallback to scipy for other dimensions or when requested
        result = sp.signal.fftconvolve(data, kernel, mode="valid")

    if normalize:
        result = result / kernel.sum()

    if axis_swap is not None:
        result = np.moveaxis(result, axis_swap[0], axis_swap[1])

    return result


def apply_directional_convolution(
    data, kernel, direction, ndim, normalize=True, use_numba=True
):
    """Apply convolution kernel for a specific direction with appropriate axis swapping.

    Parameters
    ----------
    data : np.ndarray
        Input data array
    kernel : np.ndarray
        Convolution kernel for the direction
    direction : str
        Direction ('u', 'v', or 'w')
    ndim : int
        Number of dimensions (2 or 3)
    normalize : bool, default=True
        Whether to normalize by kernel sum
    use_numba : bool, default=True
        Whether to use Numba-compiled convolution

    Returns
    -------
    np.ndarray
        Convolved result with appropriate axis swapping
    """
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
