"""JAX twin of :mod:`pybella.utils.operators.gradient`.

Corner accumulations are written as left-to-right sums in the same order as
the numpy version's sequential ``+=``, so results are bit-identical in x64.
"""

import jax.numpy as jnp

from . import finite_difference


def compute_gradient_2d(field, dx, dy):
    """Compute 2D gradient: grad(phi) = (dphi/dx, dphi/dy)."""
    grad_x = finite_difference.do_1d(field, dx, axis=0)
    grad_y = finite_difference.do_1d(field, dy, axis=1)
    return grad_x, grad_y


def compute_gradient_3d(field, dx, dy, dz):
    """Compute 3D gradient: grad(phi) = (dphi/dx, dphi/dy, dphi/dz)."""
    grad_x = finite_difference.do_1d(field, dx, axis=0)
    grad_y = finite_difference.do_1d(field, dy, axis=1)
    grad_z = finite_difference.do_1d(field, dz, axis=2)
    return grad_x, grad_y, grad_z


def _compute_grad_nodes_2d(p, dx, dy):
    """Compute 2D gradient at nodes using corner averaging."""
    bl = p[0:-1, 0:-1]  # bottom-left
    br = p[0:-1, 1:]  # bottom-right
    tl = p[1:, 0:-1]  # top-left
    tr = p[1:, 1:]  # top-right

    Dpx = -bl - br + tl + tr
    Dpy = -bl + br - tl + tr

    scale_factor = 0.5 ** (2 - 1)
    Dpx = Dpx * (scale_factor / dx)
    Dpy = Dpy * (scale_factor / dy)

    return Dpx, Dpy


def _compute_grad_nodes_3d(p, dx, dy, dz):
    """Compute 3D gradient at nodes using corner averaging."""
    c0 = p[0:-1, 0:-1, 0:-1]  # bottom-left-back
    c1 = p[0:-1, 0:-1, 1:]  # bottom-left-front
    c2 = p[0:-1, 1:, 0:-1]  # bottom-right-back
    c3 = p[0:-1, 1:, 1:]  # bottom-right-front
    c4 = p[1:, 0:-1, 0:-1]  # top-left-back
    c5 = p[1:, 0:-1, 1:]  # top-left-front
    c6 = p[1:, 1:, 0:-1]  # top-right-back
    c7 = p[1:, 1:, 1:]  # top-right-front

    Dpx = -c0 - c1 - c2 - c3 + c4 + c5 + c6 + c7
    Dpy = -c0 - c1 + c2 + c3 - c4 - c5 + c6 + c7
    Dpz = -c0 + c1 - c2 + c3 - c4 + c5 - c6 + c7

    scale_factor = 0.5 ** (3 - 1)
    Dpx = Dpx * (scale_factor / dx)
    Dpy = Dpy * (scale_factor / dy)
    Dpz = Dpz * (scale_factor / dz)

    return Dpx, Dpy, Dpz


def compute_gradient_nodes_2d(p, dx, dy):
    """Compute gradient at nodes using finite difference averaging (2D)."""
    return _compute_grad_nodes_2d(p, dx, dy)


def compute_gradient_nodes_3d(p, dx, dy, dz):
    """Compute gradient at nodes using finite difference averaging (3D)."""
    return _compute_grad_nodes_3d(p, dx, dy, dz)


def compute_at_nodes(p, ndim, dxy):
    """
    Compute gradient at nodes using finite difference averaging.

    Always returns 3 components; the z-component is zero for 2D.
    """
    dx, dy, dz = dxy

    if ndim == 2:
        grad_x, grad_y = compute_gradient_nodes_2d(p, dx, dy)
        grad_z = jnp.zeros_like(grad_x)
        return grad_x, grad_y, grad_z
    elif ndim == 3:
        return compute_gradient_nodes_3d(p, dx, dy, dz)
    else:
        raise ValueError(f"Unsupported dimension: {ndim}")
