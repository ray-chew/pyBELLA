"""Ullrich et al. (2014, 2016) analytic balanced baroclinic-wave base state.

The steady, zonally symmetric, hydrostatically balanced background used by
the Hughes & Jablonowski (2023) mountain-induced baroclinic wave test case
(GMD 16, 6805-6831; docs/gmd-16-6805-2023.pdf, Appendix B). This module
provides ONLY the DRY variant (specific humidity q_v = 0, so the virtual
temperature T_v equals the temperature T).

pyBELLA is a HEIGHT-based (Exner/theta) solver, which is the natural home
for this state: Ullrich give T, u, p and rho as closed forms of latitude
phi and geometric height z directly (Eqs. B1-B5), so — unlike the
pressure-based cores the paper targets — NO numerical root-finding of z(p)
is needed. The functions here are broadcastable over numpy arrays of
(phi, z).

The state is a westerly midlatitude jet in gradient-wind + thermal-wind
balance about a pole-to-equator temperature contrast (T_E = 310 K at the
equator, T_P = 240 K at the poles); the jet peaks near 45 deg / 10 km at
~28 m/s (paper Fig. 1a). It is the baroclinically unstable atmosphere the
two mid-latitude ridges (Eq. 1) then destabilise.

Constants are Table B1 with the small-planet factor X = 1 (full-size
Earth). ``a``, ``Omega`` carry the 1/X, X scalings so a reduced-size
variant only needs a different X.
"""

import numpy as np

# --- Table B1 physical constants (X = 1: full-size, Earth rotation) -------
X = 1.0
A_EARTH = 6.37122e6 / X  # scaled planet radius [m]
OMEGA = 2.0 * np.pi / 86164.0 * X  # scaled rotation rate [s^-1]
G = 9.80616  # gravity [m s^-2]
R_D = 287.0  # dry gas constant [J kg^-1 K^-1]
P0 = 1.0e5  # reference pressure [Pa]
B = 2.0  # jet half-width parameter
K = 3.0  # power in the temperature/jet latitudinal structure
T_E = 310.0  # equatorial surface temperature [K]
T_P = 240.0  # polar surface temperature [K]
GAMMA_LAPSE = 0.005  # vertical lapse rate [K m^-1]
CP = 1004.5  # specific heat at constant pressure [J kg^-1 K^-1]

T0 = 0.5 * (T_E + T_P)  # 275 K, the arithmetic mean surface temperature


def _height_var(z):
    """Dimensionless height A = z g / (b R_d T0) used throughout Eqs. B."""
    return z * G / (B * R_D * T0)


def _I_T(phi):
    """Latitudinal temperature-gradient structure I_T(phi) (Appendix B)."""
    c = np.cos(phi)
    return c**K - (K / (K + 2.0)) * c ** (K + 2.0)


def _tau1(z):
    a2 = _height_var(z) ** 2
    return (1.0 / T0) * np.exp(GAMMA_LAPSE * z / T0) + ((T0 - T_P) / (T0 * T_P)) * (
        1.0 - 2.0 * a2
    ) * np.exp(-a2)


def _tau2(z):
    a2 = _height_var(z) ** 2
    return (
        ((K + 2.0) / 2.0) * ((T_E - T_P) / (T_E * T_P)) * (1.0 - 2.0 * a2) * np.exp(-a2)
    )


def _tau_int1(z):
    a2 = _height_var(z) ** 2
    return (1.0 / GAMMA_LAPSE) * (np.exp(GAMMA_LAPSE * z / T0) - 1.0) + z * (
        (T0 - T_P) / (T0 * T_P)
    ) * np.exp(-a2)


def _tau_int2(z):
    a2 = _height_var(z) ** 2
    return ((K + 2.0) / 2.0) * ((T_E - T_P) / (T_E * T_P)) * z * np.exp(-a2)


def temperature(phi, z):
    """Temperature T(phi, z) [K] — Eq. B1 with q_v = 0 (T == T_v)."""
    return 1.0 / (_tau1(z) - _tau2(z) * _I_T(phi))


def pressure(phi, z):
    """Pressure p(phi, z) [Pa] — Eq. B4."""
    return P0 * np.exp(-(G / R_D) * (_tau_int1(z) - _tau_int2(z) * _I_T(phi)))


def density(phi, z):
    """Density rho(phi, z) [kg m^-3] — ideal gas, Eq. B5."""
    return pressure(phi, z) / (R_D * temperature(phi, z))


def zonal_wind(phi, z):
    """Zonal (eastward) wind u(phi, z) [m s^-1] — gradient-wind balance,
    Eq. B2 built on the auxiliary U(phi, z). Meridional wind v == 0."""
    c = np.cos(phi)
    U = (
        (G * K / A_EARTH)
        * _tau_int2(z)
        * (c ** (K - 1.0) - c ** (K + 1.0))
        * temperature(phi, z)
    )
    return -OMEGA * A_EARTH * c + np.sqrt((OMEGA * A_EARTH * c) ** 2 + A_EARTH * c * U)


# --- H&J (2023) topography: two midlatitude ridges (Sect. 2.2, Table 1) ---
# The ridges REPLACE the usual wind/temperature perturbation as the
# baroclinic-wave trigger. Gaussian in longitude, elongated (6th-power) in
# latitude, in the northern midlatitudes.
H0 = 2.0e3  # peak mountain height [m]
PHI_RIDGE = (np.pi / 4.0, np.pi / 4.0)  # ridge centre latitudes [rad]
LAM_RIDGE = (72.0 * np.pi / 180.0, 140.0 * np.pi / 180.0)  # centre longitudes [rad]
PHI_BAR = 40.0 * np.pi / 180.0  # nominal latitudinal width [rad]
LAM_BAR = 7.0 * np.pi / 180.0  # nominal longitudinal width [rad]
# scale parameters transforming the nominal widths into the Gaussian /
# 6th-power forms so z_s drops to 10% of the peak at +- half-width
_D_LAT = (PHI_BAR / 2.0) * (-np.log(0.1)) ** (-1.0 / 6.0)
_C_LON = (LAM_BAR / 2.0) * (-np.log(0.1)) ** (-1.0 / 2.0)


def _delta_lon(lam, lam_n):
    """Signed periodic longitude offset lam - lam_n wrapped to [-pi, pi).

    Its square is the paper's modified-longitude l_n(lam)^2, and its
    derivative w.r.t. lam is +-1 (Eq. 5), so it gives both the Gaussian and
    its analytic gradient with the seam handled."""
    return np.mod(lam - lam_n + np.pi, 2.0 * np.pi) - np.pi


def _ridge_terms(lam, phi):
    """The two per-ridge exponentials exp[-((dphi/d)^6 + (dlam/c)^2)] and
    their offsets, shared by the height and both gradients."""
    out = []
    for phi_n, lam_n in zip(PHI_RIDGE, LAM_RIDGE):
        dphi = phi - phi_n
        dlam = _delta_lon(lam, lam_n)
        e = np.exp(-(((dphi / _D_LAT) ** 6) + (dlam / _C_LON) ** 2))
        out.append((e, dphi, dlam))
    return out


def orography(lam, phi):
    """Surface height z_s(lam, phi) [m] — the two ridges, Eq. (1)."""
    total = 0.0
    for e, _, _ in _ridge_terms(lam, phi):
        total = total + e
    return H0 * total


def orography_grad_lam(lam, phi):
    """dz_s/dlam [m rad^-1] (analytic; SphericalTerrainMap requires it)."""
    total = 0.0
    for e, _, dlam in _ridge_terms(lam, phi):
        total = total + e * (-2.0 * dlam / _C_LON**2)
    return H0 * total


def orography_grad_phi(lam, phi):
    """dz_s/dphi [m rad^-1] (analytic)."""
    total = 0.0
    for e, dphi, _ in _ridge_terms(lam, phi):
        total = total + e * (-6.0 * (dphi / _D_LAT) ** 5 / _D_LAT)
    return H0 * total
