"""JAX twin of :mod:`pybella.flow_solver.numerics.diffusion`.

The jitted core takes the conservative fields and returns the diffused
ones; the wrapper keeps the numpy ``apply`` contract (mutates ``mem.sol``,
then refills ghost cells through the canonical numpy boundary code).
"""

import functools

import numpy as np
import jax
import jax.numpy as jnp

from ...flow_solver.utils.boundary import cell_boundary as bdry_c


def _laplacian_2d(f, dx, dy):
    interior = (f[2:, 1:-1] - 2.0 * f[1:-1, 1:-1] + f[:-2, 1:-1]) / (dx * dx) + (
        f[1:-1, 2:] - 2.0 * f[1:-1, 1:-1] + f[1:-1, :-2]
    ) / (dy * dy)
    return jnp.zeros_like(f).at[1:-1, 1:-1].set(interior)


def _laplacian_3d(f, dx, dy, dz):
    interior = (
        (f[2:, 1:-1, 1:-1] - 2.0 * f[1:-1, 1:-1, 1:-1] + f[:-2, 1:-1, 1:-1]) / (dx * dx)
        + (f[1:-1, 2:, 1:-1] - 2.0 * f[1:-1, 1:-1, 1:-1] + f[1:-1, :-2, 1:-1])
        / (dy * dy)
        + (f[1:-1, 1:-1, 2:] - 2.0 * f[1:-1, 1:-1, 1:-1] + f[1:-1, 1:-1, :-2])
        / (dz * dz)
    )
    return jnp.zeros_like(f).at[1:-1, 1:-1, 1:-1].set(interior)


@functools.partial(jax.jit, static_argnames=("ndim",))
def _diffuse(rho, rhou, rhov, rhow, rhoY, theta_bar, dtK, dxyz, ndim):
    dx, dy, dz = dxyz
    lap = (
        (lambda f: _laplacian_2d(f, dx, dy))
        if ndim == 2
        else (lambda f: _laplacian_3d(f, dx, dy, dz))
    )

    u = rhou / rho
    v = rhov / rho
    theta = rhoY / rho

    u = u + dtK * lap(u)
    v = v + dtK * lap(v)
    theta = theta + dtK * lap(theta - theta_bar)

    if ndim == 3:
        w = rhow / rho
        w = w + dtK * lap(w)

    new_rho = rhoY / theta
    new_rhou = new_rho * u
    new_rhov = new_rho * v
    new_rhow = new_rho * w if ndim == 3 else rhow
    return new_rho, new_rhou, new_rhov, new_rhow


def apply(mem, ud, dt):
    """One explicit diffusion step on velocity and theta'; updates mem.sol."""
    sol, elem = mem.sol, mem.elem
    theta_bar = 1.0 / mem.npf.HydroState.get_S0c(elem)

    rho, rhou, rhov, rhow = _diffuse(
        jnp.asarray(sol.rho),
        jnp.asarray(sol.rhou),
        jnp.asarray(sol.rhov),
        jnp.asarray(sol.rhow),
        jnp.asarray(sol.rhoY),
        jnp.asarray(theta_bar),
        dt * ud.diffusion_coeff,
        (elem.dx, elem.dy, elem.dz),
        elem.ndim,
    )
    sol.rho[...] = np.asarray(rho)
    sol.rhou[...] = np.asarray(rhou)
    sol.rhov[...] = np.asarray(rhov)
    if elem.ndim == 3:
        sol.rhow[...] = np.asarray(rhow)

    bdry_c.set_ghost_cells(mem, ud)
