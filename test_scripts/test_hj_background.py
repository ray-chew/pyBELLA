"""Ullrich (2014/2016) dry balanced background — analytic-consistency gates.

These pin the closed-form base state of the Hughes & Jablonowski (2023)
mountain-induced baroclinic wave (src/pybella/tests/ullrich_baroclinic.py)
against the paper's documented structure and its own internal balances,
BEFORE it is wired into a pyBELLA initial condition. Cheap, numpy-only.
"""

import numpy as np

from pybella.tests import ullrich_baroclinic as ub


def test_surface_temperature_endpoints():
    # equator/pole surface temperatures are the defining T_E / T_P (Table B1)
    assert abs(ub.temperature(0.0, 0.0) - ub.T_E) < 1e-9
    assert abs(ub.temperature(np.pi / 2, 0.0) - ub.T_P) < 1e-9


def test_surface_pressure_equator():
    assert abs(ub.pressure(0.0, 0.0) - ub.P0) < 1e-6


def test_midlatitude_jet_structure():
    # Fig. 1a: a westerly jet peaking near 45 deg latitude, ~10-12 km height,
    # at ~28 m s^-1; easterlies (u < 0) are absent at the jet core height
    phis = np.linspace(0.01, np.pi / 2 - 0.01, 2000)
    u = ub.zonal_wind(phis, 1.0e4)
    ipk = int(np.argmax(u))
    assert 25.0 < u[ipk] < 30.0, u[ipk]
    assert 40.0 < np.degrees(phis[ipk]) < 50.0, np.degrees(phis[ipk])
    # symmetric hemispheres (state is a function of cos(phi))
    assert np.allclose(ub.zonal_wind(phis, 1.0e4), ub.zonal_wind(-phis, 1.0e4))


def test_hydrostatic_balance():
    # dp/dz = -rho g to finite-difference accuracy at several (phi, z)
    for phi in (0.0, np.pi / 6, np.pi / 4, np.pi / 3):
        for z in (1.0e3, 5.0e3, 1.0e4, 2.0e4):
            dz = 1.0
            dpdz = (ub.pressure(phi, z + dz) - ub.pressure(phi, z - dz)) / (2 * dz)
            rhog = -ub.density(phi, z) * ub.G
            assert abs(dpdz - rhog) <= 1e-6 * abs(rhog), (phi, z, dpdz, rhog)


def test_gradient_wind_balance():
    # the zonal wind satisfies the gradient-wind relation it is built from:
    # (u/(a cos phi) + Omega) * ... i.e. u^2/(a cos) + 2 Omega u sin ... is the
    # meridional pressure-gradient force. Verify via the defining identity
    # a cos U = (u + Omega a cos)^2 - (Omega a cos)^2  (rearranged Eq. B2).
    for phi in (np.pi / 6, np.pi / 4, np.pi / 3):
        for z in (2.0e3, 1.0e4):
            c = np.cos(phi)
            u = ub.zonal_wind(phi, z)
            lhs = (u + ub.OMEGA * ub.A_EARTH * c) ** 2 - (
                ub.OMEGA * ub.A_EARTH * c
            ) ** 2
            U = (
                (ub.G * ub.K / ub.A_EARTH)
                * ub._tau_int2(z)
                * (c ** (ub.K - 1.0) - c ** (ub.K + 1.0))
                * ub.temperature(phi, z)
            )
            rhs = ub.A_EARTH * c * U
            assert abs(lhs - rhs) <= 1e-6 * max(1.0, abs(rhs)), (phi, z)


# ------------------------------------------------------- H&J ridges (Eq. 1)


def test_ridge_peaks_and_width():
    # each ridge peaks at h0 = 2000 m at its (lam_n, phi_n) centre
    for lam_n, phi_n in zip(ub.LAM_RIDGE, ub.PHI_RIDGE):
        assert abs(ub.orography(lam_n, phi_n) - ub.H0) < 1e-6
    # z_s falls to 10% of the peak at +- half nominal width (paper Sect. 2.2)
    lam_n, phi_n = ub.LAM_RIDGE[0], ub.PHI_RIDGE[0]
    assert abs(ub.orography(lam_n, phi_n + ub.PHI_BAR / 2) - 0.1 * ub.H0) < 1e-6
    assert abs(ub.orography(lam_n + ub.LAM_BAR / 2, phi_n) - 0.1 * ub.H0) < 1e-6
    # away from the ridges the surface is flat (~0)
    assert ub.orography(-np.pi / 2, -np.pi / 4) < 1e-6


def test_ridge_periodic_seam():
    # orography is periodic in longitude (the wrapped modified-longitude)
    lam_n, phi_n = ub.LAM_RIDGE[0], ub.PHI_RIDGE[0]
    assert np.allclose(
        ub.orography(lam_n, phi_n), ub.orography(lam_n + 2 * np.pi, phi_n)
    )


def test_ridge_analytic_gradients():
    # analytic dz_s/dlam, dz_s/dphi match centred finite differences
    rng = np.random.default_rng(0)
    for _ in range(20):
        lam = rng.uniform(0.0, 2 * np.pi)
        phi = rng.uniform(0.2, 1.2)  # near the northern ridges
        h = 1e-6
        fd_lam = (ub.orography(lam + h, phi) - ub.orography(lam - h, phi)) / (2 * h)
        fd_phi = (ub.orography(lam, phi + h) - ub.orography(lam, phi - h)) / (2 * h)
        assert abs(ub.orography_grad_lam(lam, phi) - fd_lam) <= 1e-4 * max(
            1.0, abs(fd_lam)
        )
        assert abs(ub.orography_grad_phi(lam, phi) - fd_phi) <= 1e-4 * max(
            1.0, abs(fd_phi)
        )
