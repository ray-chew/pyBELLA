import numpy as np
import numba as nb
from . import finite_difference


@nb.njit(cache=True)
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
    grad_x = finite_difference.do_1d(field, dx, axis=0)
    grad_y = finite_difference.do_1d(field, dy, axis=1)

    return grad_x, grad_y


@nb.njit(cache=True)
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
    grad_x = finite_difference.do_1d(field, dx, axis=0)
    grad_y = finite_difference.do_1d(field, dy, axis=1)
    grad_z = finite_difference.do_1d(field, dz, axis=2)

    return grad_x, grad_y, grad_z


#######
# Compute gradient at nodes
#######


@nb.njit(cache=True)
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


@nb.njit(cache=True)
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


@nb.njit(cache=True)
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


@nb.njit(cache=True)
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


def compute_at_nodes(p, ndim, dxy):
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
