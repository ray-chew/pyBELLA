"""Explicit constant-coefficient diffusion (e.g. for the Straka density current).

Applies forward-Euler diffusion with a fixed kinematic coefficient K to the
velocity components and to the potential-temperature *perturbation*
theta' = theta - theta_bar(y), once per full time step:

    u     <- u     + dt * K * lap(u)
    v     <- v     + dt * K * lap(v)
    theta <- theta + dt * K * lap(theta - theta_bar)

Thermal diffusion acts at fixed P = rho*theta (the pressure-like prognostic
rhoY is untouched), so the density is updated as rho = rhoY / theta — the
anomaly construction used by the bubble initial conditions.

Boundary behaviour follows from the existing ghost-cell conventions: WALL
ghosts mirror tangential velocity and scalars (zero-flux / free-slip) and
negate the normal momentum (impermeability), which is exactly the classic
Straka setup.

Stability: explicit Euler needs K*dt*(1/dx^2 + 1/dy^2 [+ 1/dz^2]) <= 1/2;
the Straka regression case sits two orders of magnitude below this.

Enabled per-case via ``ud.diffusion = True`` + ``ud.diffusion_coeff = K``
(non-dimensional: K_phys * t_ref / h_ref**2). Off by default — existing
cases are untouched.
"""

import numba as nb
import numpy as np

from ..utils.boundary import cell_boundary as bdry_c
from ...backends import is_jax_backend


@nb.njit(cache=True)
def _laplacian_2d(f, dx, dy):
    lap = np.zeros_like(f)
    lap[1:-1, 1:-1] = (f[2:, 1:-1] - 2.0 * f[1:-1, 1:-1] + f[:-2, 1:-1]) / (dx * dx) + (
        f[1:-1, 2:] - 2.0 * f[1:-1, 1:-1] + f[1:-1, :-2]
    ) / (dy * dy)
    return lap


@nb.njit(cache=True)
def _laplacian_3d(f, dx, dy, dz):
    lap = np.zeros_like(f)
    lap[1:-1, 1:-1, 1:-1] = (
        (f[2:, 1:-1, 1:-1] - 2.0 * f[1:-1, 1:-1, 1:-1] + f[:-2, 1:-1, 1:-1]) / (dx * dx)
        + (f[1:-1, 2:, 1:-1] - 2.0 * f[1:-1, 1:-1, 1:-1] + f[1:-1, :-2, 1:-1])
        / (dy * dy)
        + (f[1:-1, 1:-1, 2:] - 2.0 * f[1:-1, 1:-1, 1:-1] + f[1:-1, 1:-1, :-2])
        / (dz * dz)
    )
    return lap


def _laplacian(f, elem):
    if elem.ndim == 2:
        return _laplacian_2d(f, elem.dx, elem.dy)
    return _laplacian_3d(f, elem.dx, elem.dy, elem.dz)


def apply(mem, ud, dt):
    """One explicit diffusion step on velocity and theta'; updates mem.sol in place."""
    if is_jax_backend(ud):
        from ...backends.jax_ops import diffusion as jax_diffusion

        return jax_diffusion.apply(mem, ud, dt)

    K = ud.diffusion_coeff
    sol = mem.sol
    elem = mem.elem

    u = sol.rhou / sol.rho
    v = sol.rhov / sol.rho
    theta = sol.rhoY / sol.rho

    # background theta profile on cells (S0 = 1/theta_bar)
    theta_bar = 1.0 / mem.npf.HydroState.get_S0c(elem)

    u = u + dt * K * _laplacian(u, elem)
    v = v + dt * K * _laplacian(v, elem)
    theta = theta + dt * K * _laplacian(theta - theta_bar, elem)

    if elem.ndim == 3:
        w = sol.rhow / sol.rho
        w = w + dt * K * _laplacian(w, elem)

    # thermal diffusion at fixed P = rhoY: fold the new theta into rho
    sol.rho[...] = sol.rhoY / theta
    sol.rhou[...] = sol.rho * u
    sol.rhov[...] = sol.rho * v
    if elem.ndim == 3:
        sol.rhow[...] = sol.rho * w

    bdry_c.set_ghost_cells(mem, ud)
