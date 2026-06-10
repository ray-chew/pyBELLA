"""Smith (1980) linear hydrostatic mountain-wave reference.

Steady linear solution for uniform wind U and constant buoyancy frequency
N over a witch-of-Agnesi hill h(x) = h0 a^2 / (x^2 + a^2) in the
hydrostatic regime (N a / U >> 1), Boussinesq form:

    delta(x, z) = h0 a (a cos(l z) - x sin(l z)) / (x^2 + a^2)
    w(x, z)     = U d(delta)/dx
                = U h0 a [ (x^2 - a^2) sin(l z) - 2 a x cos(l z) ] / (x^2 + a^2)^2
    u'(x, z)    = -U d(delta)/dz
                = U h0 a l (a sin(l z) + x cos(l z)) / (x^2 + a^2)

with vertical wavenumber l = N / U, and the analytic wave drag

    D = (pi / 4) rho_0 N U h0^2        (per unit spanwise length).

The compressible simulation's wave amplitudes grow with height like
1/sqrt(rho0(z)); the comparator removes that anelastic factor by scaling
the simulated perturbations with sqrt(rho0(z) / rho0(0)) before comparing
against the Boussinesq reference.

The comparison layer (:func:`compare`) is a standalone oracle in the
spirit of ``baldauf_brdar_analytic.py``: it consumes the live ModelState
of an in-process run (driven by ``test_scripts/test_agnesi_analytic.py``)
and returns relative-L2 metrics plus the momentum-flux/drag diagnostics —
it catches *wrongness*, not just *change*.
"""

import numpy as np


def smith_fields(x, z, params):
    """Analytic delta, w, u' on a meshgrid (x[:, None], z[None, :]) in SI."""
    U, N, h0, a = params["U"], params["N"], params["h0"], params["a"]
    l = N / U
    x = np.asarray(x).reshape(-1, 1)
    z = np.asarray(z).reshape(1, -1)
    r2 = x**2 + a**2
    delta = h0 * a * (a * np.cos(l * z) - x * np.sin(l * z)) / r2
    w = U * h0 * a * ((x**2 - a**2) * np.sin(l * z) - 2 * a * x * np.cos(l * z)) / r2**2
    up = U * h0 * a * l * (a * np.sin(l * z) + x * np.cos(l * z)) / r2
    return delta, w, up


def analytic_drag(params, rho0_surface):
    """Smith wave drag per unit spanwise length [N/m]."""
    return 0.25 * np.pi * rho0_surface * params["N"] * params["U"] * params["h0"] ** 2


def _inner_xy(arr, ndim):
    """Inner-domain x-y slab; collapses the degenerate spanwise axis in 3D."""
    if ndim == 2:
        return arr[2:-2, 2:-2]
    return arr[2:-2, 2:-2, 0]


def sim_perturbations_SI(mem, ud):
    """Extract (x, z, u', w, rho0) in SI (y vertical; quasi-2D 3D or native 2D).

    Returns inner-domain cell fields with the spanwise (z-array) axis
    collapsed; heights are the physical cell heights from the metric.
    """
    ndim = mem.elem.ndim

    rho = _inner_xy(mem.sol.rho, ndim)
    u = _inner_xy(mem.sol.rhou, ndim) / rho * ud.u_ref
    w = _inner_xy(mem.sol.rhov, ndim) / rho * ud.u_ref

    x = mem.elem.x[2:-2] * ud.h_ref
    z = _inner_xy(mem.elem.metric.z, ndim) * ud.h_ref

    rho0 = _inner_xy(mem.npf.HydroState.rho0, ndim)

    up = u - ud.u_wind_speed * ud.u_ref
    # remove the anelastic 1/sqrt(rho0) amplitude growth for the
    # Boussinesq comparison
    fac = np.sqrt(rho0 / rho0[:, 0:1].mean())
    return x, z, up * fac, w * fac, rho0


def momentum_flux_profile(mem, ud, dx_SI):
    """Vertically resolved momentum flux M(eta) = -integral rho u' w dx [N/m].

    Computed on eta-levels (terrain-following), J-weighted densities; in
    the linear steady state M is height-constant and equals -D below the
    sponge.
    """
    ndim = mem.elem.ndim
    rho = _inner_xy(mem.sol.rho, ndim) * (ud.p_ref / (ud.R_gas * ud.T_ref))
    u = _inner_xy(mem.sol.rhou, ndim) / _inner_xy(mem.sol.rho, ndim) * ud.u_ref
    w = _inner_xy(mem.sol.rhov, ndim) / _inner_xy(mem.sol.rho, ndim) * ud.u_ref
    up = u - ud.u_wind_speed * ud.u_ref
    return -np.sum(rho * up * w, axis=0) * dx_SI


def compare(mem, ud, z_lo_SI=1000.0, z_hi_SI=9000.0):
    """Rel-L2 of w and u' vs Smith in the interior window + drag metrics."""
    params = {
        "U": ud.U0,
        "N": ud.NN,
        "h0": ud.hill_height,
        "a": ud.hill_width,
    }

    x, z, up_sim, w_sim, rho0 = sim_perturbations_SI(mem, ud)

    # analytic fields at the simulation's physical cell heights (column-
    # dependent under terrain): evaluate per cell
    U, N, h0, a = params["U"], params["N"], params["h0"], params["a"]
    l = N / U
    xx = x.reshape(-1, 1)
    r2 = xx**2 + a**2
    w_ref = (
        U
        * h0
        * a
        * ((xx**2 - a**2) * np.sin(l * z) - 2 * a * xx * np.cos(l * z))
        / r2**2
    )
    up_ref = U * h0 * a * l * (a * np.sin(l * z) + xx * np.cos(l * z)) / r2

    window = (z >= z_lo_SI) & (z <= z_hi_SI)

    # remove x-means from BOTH sides: on a periodic domain the wave drag
    # decelerates the mean flow during spin-up — a real effect, but not part
    # of Smith's infinite-domain steady wave solution (same convention as
    # the Baldauf-Brdar comparator's frozen k=0 mode)
    def demean(q):
        return q - q.mean(axis=0, keepdims=True)

    def rel_l2(sim, ref):
        s, r = demean(sim), demean(ref)
        denom = np.linalg.norm(np.where(window, r, 0.0))
        return np.linalg.norm(np.where(window, s - r, 0.0)) / denom

    metrics = {
        "w": rel_l2(w_sim, w_ref),
        "u": rel_l2(up_sim, up_ref),
    }

    # momentum flux vs analytic drag
    dx_SI = mem.elem.dx * ud.h_ref
    rho_ref_SI = ud.p_ref / (ud.R_gas * ud.T_ref)
    flux = momentum_flux_profile(mem, ud, dx_SI)
    drag = analytic_drag(params, rho0[:, 0].mean() * rho_ref_SI)

    zc_eta = mem.elem.y[2:-2] * ud.h_ref  # eta levels (flat away from hill)
    in_band = (zc_eta >= z_lo_SI) & (zc_eta <= z_hi_SI)
    flux_band = flux[in_band]
    metrics["drag_ratio"] = float(np.mean(flux_band) / drag)
    metrics["flux_constancy"] = float(np.std(flux_band) / np.abs(np.mean(flux_band)))

    return metrics, {"flux_profile": flux, "drag": drag, "z_eta": zc_eta}
