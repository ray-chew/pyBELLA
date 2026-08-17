# -*- coding: utf-8 -*-
import numpy as np
from numba import njit

from ....utils.operators import convolution
from ....utils import options as opts
from ....utils import slices
from ....backends import is_jax_backend

_MOMENTA = ("rhou", "rhov", "rhow")
_ALL_FLUX = ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX")


def zero_pole_faces(container, names, ig):
    """Zero the flux through the two pole faces.

    The pole face has ZERO area in the continuum, so no flux crosses it. The discrete face
    flux, built by convolving cell values, is a nonzero O(dphi^2) residual
    that BOTH pole-adjacent cells (at lambda and lambda + pi) would
    subtract with the same sign — a systematic double-loss mass/tracer
    leak. Zeroing it makes the pole a no-flux edge; over-pole transport is
    carried by the LONGITUDE fluxes of the polar cells (standard lat-lon
    FV). The face axis is the last array axis of the (sweep-oriented) flux
    container.
    """
    nfaces = getattr(container, names[0]).shape[-1]
    lo, hi = ig, nfaces - 1 - ig
    for name in names:
        arr = getattr(container, name)
        arr[..., lo] = 0.0
        arr[..., hi] = 0.0


def _normal_momentum(sol, metric, i):
    """Contravariant face momentum N_i . m for sweep axis i.

    Cartesian-component contraction, vertical component first. IEEE
    addition is not associative, so the order is fixed deliberately:
    leading with the vertical term reproduces the terrain formula
    m_v - G1 m_h1 - G2 m_h2 bit-exactly, and matches the device
    backend's separate jnp implementation of the same contraction.
    """
    Ni = metric.N[i]
    cv = metric.cart_v
    ch1, ch2 = metric.cart_haxes
    out = Ni[cv] * getattr(sol, _MOMENTA[cv]) + Ni[ch1] * getattr(sol, _MOMENTA[ch1])
    if ch2 is not None:
        out = out + Ni[ch2] * getattr(sol, _MOMENTA[ch2])
    return out


def recompute(mem, ud=None, **kwargs):
    """Recompute the advective fluxes at the cell interfaces.

    Parameters
    ----------
    mem : object
        Memory object containing sol and flux attributes
    ud : UserDataInit, optional
        When given and ud.backend == "jax", the directional convolution
        runs on the JAX backend.
    **kwargs
        Optional pre-computed velocity components ('u', 'v', 'w')
    """
    if ud is not None and is_jax_backend(ud):
        from ....backends.jax_ops import advection as jax_advection

        return jax_advection.recompute_advective_flux(mem, **kwargs)

    ndim = mem.sol.rho.ndim
    inner_idx = slices.get_inner_slice(ndim)
    kernels = convolution.get_flux_kernels(ndim)

    # Define the component order and corresponding flux indices
    components = ["u", "v"] if ndim == 2 else ["u", "v", "w"]
    rho_components = ["rhou", "rhov"] if ndim == 2 else ["rhou", "rhov", "rhow"]

    flux = mem.cache.get_flux_containers(mem.elem)

    # terrain metric; recompute is only called in the unflipped orientation
    # (time_update, between sweeps), so component i matches array axis i
    metric = mem.elem.metric

    for i, (comp, rho_comp) in enumerate(zip(components, rho_components)):
        # Use provided velocity or compute from momentum
        if comp in kwargs:
            rhoY_vel = kwargs[comp]
        else:
            momentum = getattr(mem.sol, rho_comp)
            if metric is not None:
                # general curvilinear mass flux rhoY * (N_i . m) / rho
                # (J * xi_i-dot * rhoY — what actually crosses a xi_i-face);
                # vertical-first contraction so the vertical sweep reduces
                # bit-exactly to the vertical-line contravariant flux and
                # the horizontals to the J-weighted one
                momentum = _normal_momentum(mem.sol, metric, i)
                rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho
            else:
                rhoY_vel = mem.sol.rhoY * momentum / mem.sol.rho

        # Apply directional convolution
        flux[i].rhoY[inner_idx] = convolution.apply_directional_convolution(
            rhoY_vel, kernels[comp], comp, ndim
        )

        # pole axis: the advecting mass flux through the zero-area pole face
        # must vanish, or the polar-cell reconstruction sees a spurious
        # through-pole Courant velocity
        if (
            ud is not None
            and metric is not None
            and not metric.vertical_line
            and ud.bdry_type[i] == opts.BdryType.POLE
        ):
            zero_pole_faces(flux[i], ("rhoY",), int(mem.elem.igs[i]))
