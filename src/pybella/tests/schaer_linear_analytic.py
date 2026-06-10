"""Linear mountain-wave reference for arbitrary periodic terrain (FFT).

Steady linear Boussinesq solution for uniform wind U and constant N over
ANY terrain h(x) periodic on the domain: per Fourier mode k of h,

    w_hat(k, z) = i k U h_hat(k) * E(k, z)

with the vertical structure (l = N / U the Scorer parameter)

    E = exp(i m z),  m = sign(k) sqrt(l^2 - k^2)   for |k| < l (radiating,
                                                    upward group velocity)
    E = exp(-q z),   q = sqrt(k^2 - l^2)           for |k| > l (evanescent)

and from 2D continuity u'_hat = -(1 / i k) dw_hat/dz. The Schär (2002)
ridge spectrum has its envelope band radiating hydrostatically and its
small-scale peaks (|k| = 2 pi / lambda > l) evanescent — above ~2 km the
true field carries no small-scale signal, which is the coordinate-quality
discriminator.

The wave drag is computed from the analytic surface fields,
D = -rho0 * integral u'(x,0) w(x,0) dx — convention-free; fed the Agnesi
profile it reproduces Smith's closed form D = (pi/4) rho0 N U h0^2
(self-test in test_scripts/test_schaer_analytic.py).

Mirrors ``agnesi_smith_analytic.compare`` (and reuses its simulation
extractors), adding the discriminator metric ``E_ss``: the small-scale
spectral fraction of simulated w at eta-levels where the true small-scale
signal is evanescent-dead.
"""

import numpy as np

from . import agnesi_smith_analytic as smith


def linear_fields(x_SI, z_SI, h_SI, U, N):
    """Analytic (w, u') at cells (x_i, z_ij) for periodic terrain h(x_i).

    x_SI: (Nx,) equispaced periodic cell centers; z_SI: (Nx, Nz) physical
    cell heights (column-dependent under terrain); h_SI: (Nx,) terrain.
    """
    nx = x_SI.size
    L = (x_SI[1] - x_SI[0]) * nx
    k_all = 2.0 * np.pi * np.fft.fftfreq(nx, d=L / nx)
    h_hat = np.fft.fft(h_SI) / nx  # h(x) = sum_k h_hat e^{i k x}

    l = N / U
    w = np.zeros_like(z_SI)
    up = np.zeros_like(z_SI)
    # mode loop (k = 0 carries no wave); negligible amplitudes skipped
    tiny = 1e-14 * np.abs(h_hat).max()
    for idx in range(1, nx):
        k = k_all[idx]
        a_w = 1j * k * U * h_hat[idx]
        if np.abs(a_w) < tiny:
            continue
        # DFT phases are indexed from the first sample: e^{i k (x - x0)}
        phase_x = np.exp(1j * k * (x_SI - x_SI[0])).reshape(-1, 1)
        if abs(k) < l:
            m = np.sign(k) * np.sqrt(l**2 - k**2)
            ez = np.exp(1j * m * z_SI)
            w_hat = a_w * ez
            up_hat = -(m / k) * w_hat
        else:
            q = np.sqrt(k**2 - l**2)
            ez = np.exp(-q * z_SI)
            w_hat = a_w * ez
            up_hat = -1j * (q / k) * w_hat
        w += np.real(w_hat * phase_x)
        up += np.real(up_hat * phase_x)
    return w, up


def analytic_drag(x_SI, h_SI, U, N, rho0_surface):
    """D = -rho0 integral u'(x,0) w(x,0) dx over the periodic domain [N/m]."""
    z0 = np.zeros((x_SI.size, 1))
    w0, up0 = linear_fields(x_SI, z0, h_SI, U, N)
    dx = x_SI[1] - x_SI[0]
    return -rho0_surface * np.sum(up0[:, 0] * w0[:, 0]) * dx


def small_scale_fraction(w_levels, dx_SI, lambda_SI, frac=0.75):
    """E_ss: spectral energy fraction of w at |k| >= frac * (2 pi / lambda).

    ``w_levels``: (Nx, Nlev) — w on eta-levels (k = 0 mode excluded; the
    mean carries no wave information).
    """
    nx = w_levels.shape[0]
    k = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx_SI)
    spec = np.abs(np.fft.fft(w_levels, axis=0)) ** 2
    spec[0] = 0.0
    k_cut = frac * 2.0 * np.pi / lambda_SI
    total = spec.sum()
    if total == 0.0:
        return 0.0
    return float(spec[np.abs(k) >= k_cut].sum() / total)


def compare(mem, ud, z_lo_SI=1000.0, z_hi_SI=8000.0, ss_lo_SI=4000.0, ss_hi_SI=9000.0):
    """Metrics dict in the shape of ``agnesi_smith_analytic.compare`` plus
    the coordinate-quality discriminator ``E_ss``."""
    U = ud.U0
    N = ud.NN

    x, z, up_sim, w_sim, rho0 = smith.sim_perturbations_SI(mem, ud)
    h_SI = ud.orography(x / ud.h_ref, 0.0) * ud.h_ref

    w_ref, up_ref = linear_fields(x, z, h_SI, U, N)

    window = (z >= z_lo_SI) & (z <= z_hi_SI)

    # remove x-means from both sides: periodic-domain wave drag decelerates
    # the mean flow during spin-up — real, but not part of the steady linear
    # wave (the reference k = 0 mode is zero by construction)
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

    # momentum flux vs the spectral drag
    dx_SI = mem.elem.dx * ud.h_ref
    rho_ref_SI = ud.p_ref / (ud.R_gas * ud.T_ref)
    flux = smith.momentum_flux_profile(mem, ud, dx_SI)
    drag = analytic_drag(x, h_SI, U, N, rho0[:, 0].mean() * rho_ref_SI)

    zc_eta = mem.elem.y[2:-2] * ud.h_ref  # eta levels (flat away from hill)
    in_band = (zc_eta >= z_lo_SI) & (zc_eta <= z_hi_SI)
    flux_band = flux[in_band]
    metrics["drag_ratio"] = float(np.mean(flux_band) / drag)
    metrics["flux_constancy"] = float(np.std(flux_band) / np.abs(np.mean(flux_band)))

    # discriminator: small-scale w energy at eta-levels where the true
    # small-scale response is evanescent-dead
    ss_band = (zc_eta >= ss_lo_SI) & (zc_eta <= ss_hi_SI)
    metrics["E_ss"] = small_scale_fraction(
        demean(w_sim)[:, ss_band], dx_SI, ud.ridge_wavelength
    )

    return metrics, {"flux_profile": flux, "drag": drag, "z_eta": zc_eta}
