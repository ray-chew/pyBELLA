"""Self-tests of the strict linear mountain-wave oracle (no simulation).

``tests/linear_mountain_oracle.py`` must be shown correct before it judges
the model. Independent checks:

1. isothermal background: the Exner-scaled system is constant-coefficient,
   so the numerical per-mode ODE solve must reproduce w ~ e^{lambda z}
   (exact, incl. the radiation condition) — integrator + radiation check;
2. Boussinesq limit: scaling g -> s g, T_ref -> s^2 T_ref sends H_rho,
   g / N^2 and c to infinity at fixed N; the compressible solution
   (anelastically rescaled) must converge to the closed-form Boussinesq one
   like 1 / s;
3. Smith (1980): hydrostatic Boussinesq on a large periodic domain
   reproduces Smith's closed-form Agnesi fields and drag (pi/4) rho N U h0^2;
4. Eliassen-Palm: the momentum flux rho u'w' of the compressible solution is
   height-independent (a conservation law the ODE does not impose);
5. PDE residual: the reconstructed (u', w', theta', pi') satisfy the four
   linearised compressible equations evaluated by spectral x- and
   4th-order z-differences on a fine grid (H and NH);
6. numerical convergence in the RK4 step and insensitivity to the
   radiation height z_R (the residual WKB ambiguity of the top condition).
"""

import numpy as np
import pytest

from pybella.tests import agnesi_smith_analytic as smith
from pybella.tests import linear_mountain_oracle as lo

U, N, H0, A = 10.0, 0.01, 100.0, 10e3


def _agnesi(L=200e3, nx=128, a=A):
    x = (np.arange(nx) + 0.5) * L / nx - L / 2
    return x, H0 * a**2 / (x**2 + a**2)


def _levels(nx, zmax=10e3, nz=41):
    zl = np.linspace(0.0, zmax, nz)
    return zl, np.broadcast_to(zl, (nx, nz))


def _rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


@pytest.mark.parametrize("hydrostatic", [False, True])
def test_isothermal_exact(hydrostatic):
    bg = lo.Background(isothermal=True)
    ks = np.array([3e-5, 2e-4, 1e-3, 5e-3])
    z = np.linspace(0.0, 20e3, 41)
    Y = lo.mode_structures_compressible(bg, ks, U, z, hydrostatic)
    lam, _ = lo._radiating_eig(bg, 0.0, ks, U, hydrostatic)
    for i in range(ks.size):
        ex = np.exp(lam[i] * z)
        assert np.abs(Y[0, i] - ex).max() < 1e-7 * np.abs(ex).max()


@pytest.mark.parametrize("hydrostatic", [False, True])
def test_boussinesq_limit(hydrostatic):
    x, h = _agnesi()
    zl, z = _levels(x.size)
    errs = []
    for s in (1.0, 10.0, 100.0):
        bg = lo.Background(N_const=N, T_ref=300.0 * s**2, g=9.81 * s)
        fc = lo.linear_fields(x, z, h, U, bg, "compressible", hydrostatic)
        fb = lo.linear_fields(x, z, h, U, bg, "boussinesq", hydrostatic)
        sc = np.sqrt(bg.rho(z) / bg.rho(0.0))
        errs.append(_rel(fc["w"] * sc, fb["w"]))
    # measured 2.3e-2, 8.1e-4, 7.8e-5
    assert errs[0] < 0.05
    assert errs[1] < 0.1 * errs[0] and errs[2] < 0.15 * errs[1]


def test_smith_limit_fields_and_drag():
    """Hydrostatic Boussinesq -> Smith as the periodic domain grows."""
    zl, _ = _levels(1, zmax=9e3, nz=19)
    bg = lo.Background(N_const=N)
    errs = []
    for L, nx in ((200e3, 256), (800e3, 1024), (3200e3, 4096)):
        x, h = _agnesi(L, nx)
        z = np.broadcast_to(zl, (nx, zl.size))
        f = lo.linear_fields(x, z, h, U, bg, "boussinesq", hydrostatic=True)
        _, w_s, u_s = smith.smith_fields(x, zl, {"U": U, "N": N, "h0": H0, "a": A})
        near = np.abs(x) <= 50e3
        errs.append(_rel(lo.demean(f["w"])[near], lo.demean(w_s)[near]))
    # periodic images of the 1/x^2 tails: error falls like 1/L
    assert errs[-1] < 0.01, errs
    assert errs[2] < 0.5 * errs[1] < 0.25 * errs[0], errs
    x, h = _agnesi(3200e3, 4096)
    D = lo.drag_profile(x, [0.0], h, U, bg, regime="boussinesq", hydrostatic=True)
    assert abs(D[0] / lo.smith_drag(U, N, H0, bg.rho(0.0)) - 1.0) < 0.01


@pytest.mark.parametrize("hydrostatic", [False, True])
def test_eliassen_palm_flux_constant(hydrostatic):
    x, h = _agnesi()
    bg = lo.Background(N_const=N)
    D = lo.drag_profile(
        x,
        np.linspace(0, 20e3, 21),
        h,
        U,
        bg,
        regime="compressible",
        hydrostatic=hydrostatic,
    )
    assert np.std(D) / abs(np.mean(D)) < 1e-6
    # same-domain hydrostatic Boussinesq drag (0.971 of Smith on +-100 km —
    # periodic images): compressible-H agrees to the 1/(4H^2) level
    # (measured ratio 0.9995 with the WKB radiation condition)
    if hydrostatic:
        Db = lo.drag_profile(x, [0.0], h, U, bg, regime="boussinesq", hydrostatic=True)
        assert abs(D[0] / Db[0] - 1.0) < 0.005


def _d4z(q, dz):
    """4th-order central z-derivative on interior rows."""
    return (q[:, :-4] - 8 * q[:, 1:-3] + 8 * q[:, 3:-1] - q[:, 4:]) / (12 * dz)


@pytest.mark.parametrize("hydrostatic", [False, True])
@pytest.mark.parametrize("a", [10e3, 1e3])
def test_pde_residual(hydrostatic, a):
    L, nx = 200e3, 512
    x, h = _agnesi(L, nx, a=a)
    dz = 10.0
    zl = np.arange(0.0, 8000.0 + dz / 2, dz)
    z = np.broadcast_to(zl, (nx, zl.size))
    bg = lo.Background(N_const=N)
    f = lo.linear_fields(x, z, h, U, bg, "compressible", hydrostatic)
    k = 2 * np.pi * np.fft.fftfreq(nx, d=L / nx)

    def dx_(q):
        return np.real(np.fft.ifft(1j * k[:, None] * np.fft.fft(q, axis=0), axis=0))

    c = slice(2, -2)
    zc = zl[c]
    th, pib = bg.theta(zc), bg.pi(zc)
    u, w, tp, pp = (f[q] for q in ("u", "w", "theta", "pi"))
    dnh = 0.0 if hydrostatic else 1.0
    eqs = {
        "x-mom": [U * dx_(u)[:, c], bg.cp * th * dx_(pp)[:, c]],
        "z-mom": [
            dnh * U * dx_(w)[:, c],
            bg.cp * th * _d4z(pp, dz),
            -bg.g * tp[:, c] / th,
        ],
        "theta": [U * dx_(tp)[:, c], w[:, c] * bg.dtheta_dz(zc)],
        "exner": [
            U * dx_(pp)[:, c],
            w[:, c] * bg.dpi_dz(zc),
            (bg.gamma - 1) * pib * dx_(u)[:, c],
            (bg.gamma - 1) * pib * _d4z(w, dz),
        ],
    }
    for name, terms in eqs.items():
        res = np.linalg.norm(sum(terms))
        scale = max(np.linalg.norm(t) for t in terms)
        assert res < 1e-5 * scale, f"{name}: {res / scale:.2e}"


@pytest.mark.parametrize("hydrostatic", [False, True])
@pytest.mark.parametrize("a", [10e3, 1e3])
def test_convergence_and_radiation_height(hydrostatic, a):
    """RK4 step convergence and z_R-insensitivity (default z_R = 24 km).

    Measured (w rel-L2 over 0-12 km): dz 5 m vs 2.5 m ~4e-11; z_R in
    12-20 km ~3e-4, except the narrow NH hill (~1.5e-2: modes with
    k ~ l have vertical wavelengths of tens of km and turning points,
    where no WKB radiation condition is unambiguous).
    """
    x, h = _agnesi(nx=256, a=a)
    zl, z = _levels(x.size, zmax=12e3, nz=49)
    bg = lo.Background(N_const=N)
    ref = lo.linear_fields(x, z, h, U, bg, hydrostatic=hydrostatic, dz=2.5)
    f5 = lo.linear_fields(x, z, h, U, bg, hydrostatic=hydrostatic, dz=5.0)
    assert _rel(f5["w"], ref["w"]) < 1e-8
    tol = 3e-2 if (a < 5e3 and not hydrostatic) else 1e-3
    for zR in (12e3, 16e3, 20e3):
        fR = lo.linear_fields(x, z, h, U, bg, hydrostatic=hydrostatic, z_R=zR)
        assert _rel(fR["w"], ref["w"]) < tol
