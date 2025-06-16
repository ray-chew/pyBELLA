import numpy as np
import scipy as sp
from numba import njit
from functools import lru_cache

@lru_cache(maxsize=2)
def create_convolution_kernels(ndim):
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
        return {
            'u': kernel_u,
            'v': kernel_u.T
        }
    elif ndim == 3:
        kernel_u = np.array([
            [[1, 2, 1], [2, 4, 2], [1, 2, 1]], 
            [[1, 2, 1], [2, 4, 2], [1, 2, 1]]
        ])
        return {
            'u': kernel_u,
            'v': np.swapaxes(kernel_u, 1, 0),
            'w': np.swapaxes(kernel_u, 2, 0)
        }
    else:
        raise ValueError(f"Unsupported dimension: {ndim}")


@njit
def _numba_convolve_2d(data, kernel):
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


@njit
def _numba_convolve_3d(data, kernel):
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
                            result[i, j, k] += data[i + ki, j + kj, k + kk] * kernel[ki, kj, kk]
    
    return result


def apply_convolution_kernel(data, kernel, normalize=True, axis_swap=None, use_numba=True):
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
            result = _numba_convolve_2d(data, kernel)
        else:  # 3D
            result = _numba_convolve_3d(data, kernel)
    else:
        # Fallback to scipy for other dimensions or when requested
        result = sp.signal.fftconvolve(data, kernel, mode="valid")
    
    if normalize:
        result = result / kernel.sum()
    
    if axis_swap is not None:
        result = np.moveaxis(result, axis_swap[0], axis_swap[1])
    
    return result


# Configuration for directional convolutions
DIRECTION_CONFIG = {
    2: {
        'u': {'axis_swap': (0, -1)},
        'v': {'axis_swap': None}
    },
    3: {
        'u': {'axis_swap': (0, -1)},
        'v': {'axis_swap': (-1, 0)},
        'w': {'axis_swap': None}
    }
}

def apply_directional_convolution(data, kernel, direction, ndim, normalize=True, use_numba=True):
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
        data, kernel, 
        normalize=normalize, 
        axis_swap=config['axis_swap'], 
        use_numba=use_numba
    )


@njit(cache=True)
def compute_divergence_2d(u_field, v_field, dx, dy):
    """
    Compute 2D divergence: ∇·F = ∂u/∂x + ∂v/∂y
    
    Parameters
    ----------
    u_field : np.ndarray
        Field component in x-direction
    v_field : np.ndarray  
        Field component in y-direction
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
        
    Returns
    -------
    np.ndarray
        Divergence field averaged to cell centers
    """
    # X-direction: ∂u/∂x
    div_x = finite_difference_1d(u_field, dx, axis=0)
    # Average to y-cell centers
    div_x = 0.5 * (div_x[:, :-1] + div_x[:, 1:])
    
    # Y-direction: ∂v/∂y  
    div_y = finite_difference_1d(v_field, dy, axis=1)
    # Average to x-cell centers
    div_y = 0.5 * (div_y[:-1, :] + div_y[1:, :])
    
    return div_x + div_y


@njit(cache=True)
def compute_divergence_3d(u_field, v_field, w_field, dx, dy, dz):
    """
    Compute 3D divergence: ∇·F = ∂u/∂x + ∂v/∂y + ∂w/∂z
    
    Parameters
    ----------
    u_field : np.ndarray
        Field component in x-direction
    v_field : np.ndarray
        Field component in y-direction  
    w_field : np.ndarray
        Field component in z-direction
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
    dz : float
        Grid spacing in z-direction
        
    Returns
    -------
    tuple
        (div_x, div_y, div_z) - Individual divergence components
    """
    # X-direction: ∂u/∂x
    div_x = finite_difference_1d(u_field, dx, axis=0)
    # Average to y-cell centers, then to z-faces
    div_x = 0.5 * (div_x[:, :-1, :] + div_x[:, 1:, :])
    div_x = -0.5 * (div_x[:, :, :-1] + div_x[:, :, 1:])  # Note: negative from original
    
    # Y-direction: ∂v/∂y
    div_y = finite_difference_1d(v_field, dy, axis=1)
    # Average to x-cell centers, then to z-faces
    div_y = 0.5 * (div_y[:-1, :, :] + div_y[1:, :, :])
    div_y = 0.5 * (div_y[:, :, :-1] + div_y[:, :, 1:])
    
    # Z-direction: ∂w/∂z
    div_z = finite_difference_1d(w_field, dz, axis=2)
    # Average to cell centers
    div_z = 0.5 * (div_z[:-1, :, :] + div_z[1:, :, :])
    div_z = 0.5 * (div_z[:, :-1, :] + div_z[:, 1:, :])
    
    return div_x, div_y, div_z


@njit(cache=True)
def compute_divergence_3d_total(u_field, v_field, w_field, dx, dy, dz):
    """
    Compute total 3D divergence: ∇·F = ∂u/∂x + ∂v/∂y + ∂w/∂z
    
    Parameters
    ----------
    u_field : np.ndarray
        Field component in x-direction
    v_field : np.ndarray
        Field component in y-direction  
    w_field : np.ndarray
        Field component in z-direction
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
    dz : float
        Grid spacing in z-direction
        
    Returns
    -------
    np.ndarray
        Total divergence field
    """
    div_x, div_y, div_z = compute_divergence_3d(u_field, v_field, w_field, dx, dy, dz)
    return div_x + div_y + div_z


@njit(cache=True)
def finite_difference_1d(field, spacing, axis=0):
    """
    Compute 1D finite difference along specified axis.
    
    Parameters
    ----------
    field : np.ndarray
        Input field
    spacing : float
        Grid spacing
    axis : int, default=0
        Axis along which to compute difference
        
    Returns
    -------
    np.ndarray
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
    
# @njit
def average_to_centers_2d(field, axis):
    """
    Average field values to cell centers along specified axis.
    
    Parameters
    ----------
    field : np.ndarray
        Input field (2D)
    axis : int
        Axis along which to average (0 or 1)
        
    Returns
    -------
    np.ndarray
        Averaged field
    """
    if axis == 0:
        return 0.5 * (field[:-1, :] + field[1:, :])
    elif axis == 1:
        return 0.5 * (field[:, :-1] + field[:, 1:])
    else:
        raise ValueError("axis must be 0 or 1 for 2D arrays")


# @njit
def average_to_centers_3d(field, axis):
    """
    Average field values to cell centers along specified axis.
    
    Parameters
    ----------
    field : np.ndarray
        Input field (3D)
    axis : int
        Axis along which to average (0, 1, or 2)
        
    Returns
    -------
    np.ndarray
        Averaged field
    """
    if axis == 0:
        return 0.5 * (field[:-1, :, :] + field[1:, :, :])
    elif axis == 1:
        return 0.5 * (field[:, :-1, :] + field[:, 1:, :])
    elif axis == 2:
        return 0.5 * (field[:, :, :-1] + field[:, :, 1:])
    else:
        raise ValueError("axis must be 0, 1, or 2 for 3D arrays")


@njit(cache=True)
def compute_gradient_2d(field, dx, dy):
    """
    Compute 2D gradient: ∇φ = (∂φ/∂x, ∂φ/∂y)
    
    Parameters
    ----------
    field : np.ndarray
        Scalar field
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
        
    Returns
    -------
    tuple
        (grad_x, grad_y) - Gradient components
    """
    grad_x = finite_difference_1d(field, dx, axis=0)
    grad_y = finite_difference_1d(field, dy, axis=1)
    
    return grad_x, grad_y


@njit(cache=True)
def compute_gradient_3d(field, dx, dy, dz):
    """
    Compute 3D gradient: ∇φ = (∂φ/∂x, ∂φ/∂y, ∂φ/∂z)
    
    Parameters
    ----------
    field : np.ndarray
        Scalar field
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
    dz : float
        Grid spacing in z-direction
        
    Returns
    -------
    tuple
        (grad_x, grad_y, grad_z) - Gradient components
    """
    grad_x = finite_difference_1d(field, dx, axis=0)
    grad_y = finite_difference_1d(field, dy, axis=1)
    grad_z = finite_difference_1d(field, dz, axis=2)
    
    return grad_x, grad_y, grad_z


#######
# Compute gradient at nodes
#######

@njit(cache=True)
def _compute_grad_nodes_2d(p, dx, dy):
    """Compute 2D gradient at nodes using corner averaging."""
    # Pre-computed signs for each corner
    signs_x = np.array([-1.0, -1.0, +1.0, +1.0])
    signs_y = np.array([-1.0, +1.0, -1.0, +1.0])
    
    # Initialize gradient components
    Dpx = np.zeros((p.shape[0] - 1, p.shape[1] - 1))
    Dpy = np.zeros((p.shape[0] - 1, p.shape[1] - 1))
    
    # Corner contributions
    # Bottom-left
    Dpx += signs_x[0] * p[0:-1, 0:-1]
    Dpy += signs_y[0] * p[0:-1, 0:-1]
    
    # Bottom-right
    Dpx += signs_x[1] * p[0:-1, 1:]
    Dpy += signs_y[1] * p[0:-1, 1:]
    
    # Top-left
    Dpx += signs_x[2] * p[1:, 0:-1]
    Dpy += signs_y[2] * p[1:, 0:-1]
    
    # Top-right
    Dpx += signs_x[3] * p[1:, 1:]
    Dpy += signs_y[3] * p[1:, 1:]
    
    # Apply scaling factors
    scale_factor = 0.5 ** (2 - 1)  # 0.5^(ndim-1)
    Dpx *= scale_factor / dx
    Dpy *= scale_factor / dy
    
    return Dpx, Dpy

@njit(cache=True)
def _compute_grad_nodes_3d(p, dx, dy, dz):
    """Compute 3D gradient at nodes using corner averaging."""
    # Pre-computed signs for each corner
    signs_x = np.array([-1.0, -1.0, -1.0, -1.0, +1.0, +1.0, +1.0, +1.0])
    signs_y = np.array([-1.0, -1.0, +1.0, +1.0, -1.0, -1.0, +1.0, +1.0])
    signs_z = np.array([-1.0, +1.0, -1.0, +1.0, -1.0, +1.0, -1.0, +1.0])
    
    # Initialize gradient components
    Dpx = np.zeros((p.shape[0] - 1, p.shape[1] - 1, p.shape[2] - 1))
    Dpy = np.zeros((p.shape[0] - 1, p.shape[1] - 1, p.shape[2] - 1))
    Dpz = np.zeros((p.shape[0] - 1, p.shape[1] - 1, p.shape[2] - 1))
    
    # Corner contributions (8 corners for 3D)
    # Bottom-left-back
    Dpx += signs_x[0] * p[0:-1, 0:-1, 0:-1]
    Dpy += signs_y[0] * p[0:-1, 0:-1, 0:-1]
    Dpz += signs_z[0] * p[0:-1, 0:-1, 0:-1]
    
    # Bottom-left-front
    Dpx += signs_x[1] * p[0:-1, 0:-1, 1:]
    Dpy += signs_y[1] * p[0:-1, 0:-1, 1:]
    Dpz += signs_z[1] * p[0:-1, 0:-1, 1:]
    
    # Bottom-right-back
    Dpx += signs_x[2] * p[0:-1, 1:, 0:-1]
    Dpy += signs_y[2] * p[0:-1, 1:, 0:-1]
    Dpz += signs_z[2] * p[0:-1, 1:, 0:-1]
    
    # Bottom-right-front
    Dpx += signs_x[3] * p[0:-1, 1:, 1:]
    Dpy += signs_y[3] * p[0:-1, 1:, 1:]
    Dpz += signs_z[3] * p[0:-1, 1:, 1:]
    
    # Top-left-back
    Dpx += signs_x[4] * p[1:, 0:-1, 0:-1]
    Dpy += signs_y[4] * p[1:, 0:-1, 0:-1]
    Dpz += signs_z[4] * p[1:, 0:-1, 0:-1]
    
    # Top-left-front
    Dpx += signs_x[5] * p[1:, 0:-1, 1:]
    Dpy += signs_y[5] * p[1:, 0:-1, 1:]
    Dpz += signs_z[5] * p[1:, 0:-1, 1:]
    
    # Top-right-back
    Dpx += signs_x[6] * p[1:, 1:, 0:-1]
    Dpy += signs_y[6] * p[1:, 1:, 0:-1]
    Dpz += signs_z[6] * p[1:, 1:, 0:-1]
    
    # Top-right-front
    Dpx += signs_x[7] * p[1:, 1:, 1:]
    Dpy += signs_y[7] * p[1:, 1:, 1:]
    Dpz += signs_z[7] * p[1:, 1:, 1:]
    
    # Apply scaling factors
    scale_factor = 0.5 ** (3 - 1)  # 0.5^(ndim-1)
    Dpx *= scale_factor / dx
    Dpy *= scale_factor / dy
    Dpz *= scale_factor / dz
    
    return Dpx, Dpy, Dpz

@njit(cache=True)
def compute_gradient_nodes_2d(p, dx, dy):
    """
    Compute gradient at nodes using finite difference averaging.
    
    Parameters
    ----------
    p : np.ndarray
        Scalar field (2D)
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
        
    Returns
    -------
    tuple
        (grad_x, grad_y) - Gradient components at nodes
        
    Notes
    -----
    The gradient is computed at nodes by averaging contributions from
    all neighboring cells. Output arrays have shape (nx-1, ny-1).
    """
    return _compute_grad_nodes_2d(p, dx, dy)

@njit(cache=True)
def compute_gradient_nodes_3d(p, dx, dy, dz):
    """
    Compute gradient at nodes using finite difference averaging.
    
    Parameters
    ----------
    p : np.ndarray
        Scalar field (3D)
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
    dz : float
        Grid spacing in z-direction
        
    Returns
    -------
    tuple
        (grad_x, grad_y, grad_z) - Gradient components at nodes
        
    Notes
    -----
    The gradient is computed at nodes by averaging contributions from
    all neighboring cells. Output arrays have shape (nx-1, ny-1, nz-1).
    """
    return _compute_grad_nodes_3d(p, dx, dy, dz)

def compute_gradient_nodes(p, ndim, dxy):
    """
    Compute gradient at nodes using finite difference averaging.
    
    Parameters
    ----------
    p : np.ndarray
        Scalar field
    ndim : int
        Number of dimensions (2 or 3)
    dxy : tuple
        Grid spacings (dx, dy, dz)
        
    Returns
    -------
    tuple
        Gradient components at nodes. For 2D: (grad_x, grad_y, zeros)
        For 3D: (grad_x, grad_y, grad_z)
        
    Notes
    -----
    This is the main interface function that dispatches to the appropriate
    dimension-specific implementation. Always returns 3 components for
    consistency, with the z-component being zero for 2D cases.
    """
    dx, dy, dz = dxy
    
    if ndim == 2:
        grad_x, grad_y = compute_gradient_nodes_2d(p, dx, dy)
        grad_z = np.zeros_like(grad_x)
        return grad_x, grad_y, grad_z
    elif ndim == 3:
        return compute_gradient_nodes_3d(p, dx, dy, dz)
    else:
        raise ValueError(f"Unsupported dimension: {ndim}")