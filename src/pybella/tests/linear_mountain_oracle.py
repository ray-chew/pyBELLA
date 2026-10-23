r"""Strict steady linear mountain-wave oracle for periodic terrain.

Steady, inviscid, adiabatic, 2D (x, z) flow of uniform wind U over periodic
terrain h(x), linearised about the *model's own* hydrostatic background, on
the model's actual periodic domain. Three regimes share one interface:

* ``"boussinesq"``      Boussinesq, constant N (closed form per Fourier mode),
* ``"compressible"``    fully compressible Euler about the model background
                        (per-mode vertical ODE, solved numerically),

each with ``hydrostatic=True/False`` (drop / keep the vertical acceleration
``D w / D t``). ``boussinesq + hydrostatic`` on an infinite domain is
Smith (1980); ``compressible + nonhydrostatic`` is what pyBELLA solves with
``is_compressible = 1, is_nonhydrostatic = 1``; ``compressible +
hydrostatic`` is the linear target of ``is_nonhydrostatic = 0``.

Background (derived from the code, not assumed)
-----------------------------------------------
``hydrostatics.integrated_state`` (terrain branch
``_integrated_state_fields``) builds, at physical height z,

    theta(z) = T_ref * S(z)^{-1},       S^{-1} = ud.stratification(z / h_ref)
    pi(z)    = 1 - (g / (c_p T_ref)) * int_0^z T_ref / theta dz'
    P = rho theta = (p_ref / R) pi^{1/(gamma-1)},   rho = P / theta.

For the mountain cases ``stratification = exp(N^2 z / g)``: the background
is **constant-N, NOT isothermal**. T = theta pi falls from 300 K at the
ground to ~125 K at 24 km, and pi -> 0 at z* = -(g/N^2) ln(1 - c_p T_ref
N^2 / g^2) ~ 37 km (the constant-N atmosphere has a finite top). The
density scale height H_rho = -(d ln rho / dz)^{-1} is ~10.9 km at the
surface and shrinks with height, so no constant-coefficient formula is
exact; :class:`Background` evaluates theta, pi, N^2 = g theta'/theta,
c^2 = gamma R theta pi and rho to quadrature accuracy (closed form for
constant N).

Linearised equations (pyBELLA's variables: u, w, theta, Exner pi)
-----------------------------------------------------------------
With primes for perturbations and c_p theta grad pi the pressure force,

    (1) U u'_x + c_p thb pi'_x                    = 0
    (2) delta_nh U w'_x + c_p thb pi'_z - g theta'/thb = 0
    (3) U theta'_x + w' thb_z                     = 0
    (4) U pi'_x + w' pib_z + (R/c_v) pib (u'_x + w'_z) = 0

((4) is D pi / D t = -(gamma - 1) pi div u, i.e. the P = rho theta
equation pyBELLA advances; delta_nh = 1 nonhydrostatic, 0 hydrostatic).
For a mode e^{ikx} and c^2 = gamma R thb pib, eliminating u' and theta'
leaves the first-order system for Y = (w_hat, pi_hat):

    w_hat' = -[ pib_z w_hat + i k (U - c^2/U) pi_hat ] / ((gamma-1) pib)
    pi_hat' = -i (delta_nh k U - N^2 / (k U)) w_hat / (c_p thb)

and u_hat = -c_p thb pi_hat / U, theta_hat = i thb_z w_hat / (k U).
In the anelastic limit (c -> inf) with w_hat = sqrt(rho(0)/rho(z)) W this is

    W'' + [ N^2/U^2 - delta_nh k^2 - 1/(4 H_rho^2) + (1/(2 H_rho))' ] W = 0,

i.e. the 1/sqrt(rho) amplitude growth, the 1/(4H^2) term and the
(1/2H)' correction (non-zero here: constant N is not isothermal) — plus
O(U^2 k^2 / c^2) acoustic corrections that the full system keeps. The
Boussinesq limit drops all three (rho, H -> const, c -> inf):
w_hat'' + (l^2 - delta_nh k^2) w_hat = 0, l = N / U.

Boundary and radiation conditions
---------------------------------
Bottom: linearised free slip at z = 0, w_hat(0) = i k U h_hat (the model
applies it on z = h(x); the difference is O(h^2), i.e. nonlinear).
Top: the infinite-atmosphere radiation condition. Boussinesq: exact,
w ~ e^{i m z}, m = sgn(kU) sqrt(l^2 - delta k^2) for propagating modes
(upward group velocity for U > 0 needs sgn m = sgn k), e^{-q z} for
evanescent ones. Compressible: the background varies (and the constant-N
atmosphere ends at z* ~ 37 km), so an exact radiation condition does not
exist; it is imposed at a radiation height z_R (default 24 km, the model
top) as the first-order Liouville-Green (WKB) condition on the
vertical-velocity equation (see :func:`_radiating_eig`) — upward group
velocity, or decay, *including* the Q^{-1/4} amplitude variation. A
frozen-coefficient eigenvector instead leaves a z_R-dependent partial
reflection worth +-3 % in the drag (measured); with the WKB condition the
surface drag varies by < 0.03 % for z_R in 8-24 km. The system is then
integrated *downward* to z = 0 (classical RK4, dz = 5 m, renormalised
every step; the stable direction for evanescent modes) and normalised to
the bottom BC. Residual ambiguity (``test_scripts/test_linear_oracle.py``):
~3e-4 rel-L2 in w for z_R in 12-20 km, except near-cutoff modes
(k ~ l, vertical wavelengths of tens of km, turning points) of a narrow
nonhydrostatic hill, ~1.5e-2.

k = 0 carries no wave (the periodic mean-flow deceleration by the drag is
a secular, not a steady linear, effect) — comparisons x-demean both sides.

Diagnostics
-----------
``drag``: D = -int rho_b u' w' dx at z (Eliassen-Palm: height-independent
in the steady linear compressible solution — a self-test). For the
Boussinesq regimes rho_b = rho(0) is used.
"""

import numpy as np

REGIMES = ("boussinesq", "compressible")


class Background:
    """Dimensional hydrostatic background as a function of physical height.

    ``strat(z_SI)`` returns theta / theta0 (theta0 = T_ref). If ``N_const``
    is given the profile is taken as exp(N^2 z / g) with closed-form Exner
    pressure; otherwise pi is computed by adaptive-free fine quadrature.
    """

    def __init__(
        self,
        R=287.4,
        gamma=1.4,
        T_ref=300.0,
        p_ref=1e5,
        g=9.81,
        N_const=None,
        strat=None,
        isothermal=False,
        z_max=40e3,
    ):
        self.R, self.gamma, self.T_ref, self.p_ref, self.g = R, gamma, T_ref, p_ref, g
        self.cp = gamma * R / (gamma - 1.0)
        self.N_const = N_const
        self.isothermal = isothermal
        self._strat = strat
        if strat is not None and N_const is None and not isothermal:
            zf = np.linspace(0.0, z_max, 400001)
            inv = 1.0 / strat(zf)
            integ = np.concatenate(
                ([0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(zf)))
            )
            self._zf, self._integ = zf, integ

    @classmethod
    def from_ud(cls, ud):
        """The model background: same constants + ud.stratification.

        Recognises the mountain-case constant-N profile (ud.NN set and
        stratification == exp(N^2 z / g)) and then uses the closed form.
        """
        h_ref = ud.h_ref
        strat = lambda z: ud.stratification(np.asarray(z) / h_ref)  # noqa: E731
        NN = getattr(ud, "NN", None)
        kw = dict(R=ud.R_gas, gamma=ud.gamm, T_ref=ud.T_ref, p_ref=ud.p_ref, g=ud.grav)
        if NN is not None:
            zt = np.linspace(0.0, 30e3, 7)
            if np.allclose(strat(zt), np.exp(NN**2 * zt / ud.grav), rtol=1e-12):
                return cls(N_const=NN, **kw)
        return cls(strat=strat, **kw)

    # --- profiles -----------------------------------------------------------

    def theta(self, z):
        z = np.asarray(z, dtype=float)
        if self.isothermal:
            return self.T_ref / self.pi(z)
        if self.N_const is not None:
            return self.T_ref * np.exp(self.N_const**2 * z / self.g)
        return self.T_ref * self._strat(z)

    def pi(self, z):
        z = np.asarray(z, dtype=float)
        if self.isothermal:
            # T = T_ref const: d pi/dz = -g pi / (c_p T_ref)
            return np.exp(-self.g * z / (self.cp * self.T_ref))
        if self.N_const is not None:
            a = self.N_const**2 / self.g
            integ = (1.0 - np.exp(-a * z)) / a
        else:
            integ = np.interp(z, self._zf, self._integ)
        return 1.0 - self.g / (self.cp * self.T_ref) * integ

    def dpi_dz(self, z):
        return -self.g / (self.cp * self.theta(z))

    def N2(self, z):
        z = np.asarray(z, dtype=float)
        if self.isothermal:
            return self.g**2 / (self.cp * self.T_ref) + 0.0 * z
        if self.N_const is not None:
            return self.N_const**2 + 0.0 * z
        dz = 1e-3
        return (
            self.g
            * (np.log(self.theta(z + dz)) - np.log(self.theta(z - dz)))
            / (2 * dz)
        )

    def dtheta_dz(self, z):
        return self.theta(z) * self.N2(z) / self.g

    def c2(self, z):
        return self.gamma * self.R * self.theta(z) * self.pi(z)

    def rho(self, z):
        pi = self.pi(z)
        return self.p_ref * pi ** (1.0 / (self.gamma - 1.0)) / (self.R * self.theta(z))


# --- per-mode vertical structure -------------------------------------------


def _coeffs(bg, z, ks, U, hydrostatic):
    """Entries of A(z) in Y' = A Y, Y = (w_hat, pi_hat); shape (nk,).

    A = [[A11, A12], [A21, 0]] (see module docstring).
    """
    pib = bg.pi(z)
    a = (bg.gamma - 1.0) * pib
    A11 = -bg.dpi_dz(z) / a
    A12 = -1j * ks * (U - bg.c2(z) / U) / a
    dnh = 0.0 if hydrostatic else 1.0
    A21 = -1j * (dnh * ks * U - bg.N2(z) / (ks * U)) / (bg.cp * bg.theta(z))
    return A11, A12, A21


def _pq(bg, z, ks, U, hydrostatic):
    """w'' + p w' + q w = 0 equivalent of Y' = A Y (pi = (w' - A11 w)/A12)."""
    d = 1.0  # [m] centred differences of smooth background coefficients
    A11, A12, A21 = _coeffs(bg, z, ks, U, hydrostatic)
    A11p, A12p, _ = _coeffs(bg, z + d, ks, U, hydrostatic)
    A11m, A12m, _ = _coeffs(bg, z - d, ks, U, hydrostatic)
    dA11 = (A11p - A11m) / (2 * d)
    lA12 = ((A12p - A12m) / (2 * d) / A12).real  # A12'/A12 is real
    p = -(A11 + lA12)
    q = (-dA11 + A11 * lA12 - A12 * A21).real
    return p, q


def _radiating_eig(bg, z, ks, U, hydrostatic):
    """Radiation (or decay) condition at height z, per k: (w'/w, (w, pi)).

    Liouville-Green (first-order WKB) condition. With w = W exp(-int p/2)
    the vertical-velocity equation is W'' + Q W = 0,
    Q = q - p'/2 - p^2/4, and the upward-radiating WKB solution
    W ~ Q^{-1/4} exp(i sgn(kU) int sqrt(Q)) gives

        w'/w = i sgn(kU) sqrt(Q) - Q'/(4Q) - p/2        (Q > 0)
        w'/w = -sqrt(-Q) - Q'/(4Q) - p/2                (Q < 0, decaying)

    (sgn(kU): upward group velocity for steady lee waves). This is exact
    for an isothermal background (Q, p const) and leaves only O(second-
    order WKB) reflection otherwise — far smaller than a frozen-coefficient
    eigenvector, which ignores the background's variation at z_R and
    leaves a z_R-dependent partial reflection of several % in the drag.
    pi follows from w' = A11 w + A12 pi.
    """
    d = 1.0
    p, q = _pq(bg, z, ks, U, hydrostatic)
    pp, qp = _pq(bg, z + d, ks, U, hydrostatic)
    pm, qm = _pq(bg, z - d, ks, U, hydrostatic)
    dp = (pp - pm) / (2 * d)
    Q = q - 0.5 * dp - 0.25 * p**2
    Qp = qp - 0.5 * dp - 0.25 * pp**2
    Qm = qm - 0.5 * dp - 0.25 * pm**2
    # Q' (p'' contributions are O(d^2)-small for smooth backgrounds)
    dQ = (Qp - Qm) / (2 * d)
    sq = np.sqrt(np.abs(Q))
    lam = np.where(Q > 0, 1j * np.sign(ks * U) * sq, -sq + 0j)
    # amplitude (Q^{-1/4}) correction only where WKB is valid; near a
    # turning point (Q ~ 0, |Q'| / |Q|^{3/2} large) keep the frozen root
    wkb_ok = np.abs(dQ) < 0.5 * np.abs(Q) ** 1.5
    lam = lam - np.where(wkb_ok, dQ / (4 * Q), 0.0) - 0.5 * p
    A11, A12, _ = _coeffs(bg, z, ks, U, hydrostatic)
    return lam, np.stack([np.ones_like(lam), (lam - A11) / A12])


def mode_structures_compressible(
    bg, ks, U, z_eval, hydrostatic=False, z_R=24e3, dz=5.0
):
    """Vertical structures (w_hat, pi_hat)(z) for unit w_hat(0), all k at once.

    Classical RK4 with fixed step dz, integrated downward from z_R (the
    radiation height, WKB radiation condition) to min(z_eval, 0),
    renormalised every step (evanescent modes grow by e^{q z_R} downward),
    then cubic-Hermite interpolated (using Y' = A Y) to z_eval.
    Returns array (2, nk, len(z_eval)). Points above z_R continue the
    local WKB solution exp(lambda (z - z_R)).
    """
    ks = np.asarray(ks, dtype=float)
    z_eval = np.asarray(z_eval, dtype=float).ravel()
    z_lo = min(0.0, float(z_eval.min()))
    n1 = int(np.ceil(z_R / dz))
    h = -z_R / n1  # < 0; z = 0 lands exactly on the grid
    n2 = int(np.ceil(z_lo / h)) if z_lo < 0 else 0
    n = n1 + n2
    zg = z_R + h * np.arange(n + 1)  # descending
    zg[n1] = 0.0

    def f(z, Y):
        A11, A12, A21 = _coeffs(bg, z, ks, U, hydrostatic)
        return np.stack([A11 * Y[0] + A12 * Y[1], A21 * Y[0]])

    lam_R, Y = _radiating_eig(bg, z_R, ks, U, hydrostatic)
    Y = Y / np.abs(Y).max(axis=0)
    Ys = np.empty((n + 1, 2, ks.size), dtype=complex)
    logs = np.zeros((n + 1, ks.size))
    Ys[0] = Y
    for j in range(n):
        z = zg[j]
        k1 = f(z, Y)
        k2 = f(z + 0.5 * h, Y + 0.5 * h * k1)
        k3 = f(z + 0.5 * h, Y + 0.5 * h * k2)
        k4 = f(z + h, Y + h * k3)
        Y = Y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        s = np.abs(Y).max(axis=0)
        Y = Y / s
        logs[j + 1] = logs[j] + np.log(s)
        Ys[j + 1] = Y
    # normalise to w_hat(0) = 1
    i0 = n1
    scale = np.exp(logs - logs[i0])[:, None, :] / Ys[i0, 0][None, None, :]
    S = Ys * scale  # (n+1, 2, nk)
    dS = np.empty_like(S)
    for j in range(n + 1):
        dS[j] = f(zg[j], S[j])

    out = np.empty((2, ks.size, z_eval.size), dtype=complex)
    above = z_eval > z_R
    zb = z_eval[~above]
    # ascending view for searchsorted
    za = zg[::-1]
    Sa, dSa = S[::-1], dS[::-1]
    j = np.clip(np.searchsorted(za, zb) - 1, 0, n - 1)
    hh = za[j + 1] - za[j]
    t = ((zb - za[j]) / hh)[None, None, :]
    h00 = 2 * t**3 - 3 * t**2 + 1
    h10 = t**3 - 2 * t**2 + t
    h01 = -2 * t**3 + 3 * t**2
    h11 = t**3 - t**2
    P0 = np.moveaxis(Sa[j], 0, -1)  # (2, nk, nb)
    P1 = np.moveaxis(Sa[j + 1], 0, -1)
    D0 = np.moveaxis(dSa[j], 0, -1)
    D1 = np.moveaxis(dSa[j + 1], 0, -1)
    hb = hh[None, None, :]
    out[:, :, ~above] = h00 * P0 + h10 * hb * D0 + h01 * P1 + h11 * hb * D1
    if np.any(above):
        dzu = z_eval[above] - z_R
        out[:, :, above] = S[0][:, :, None] * np.exp(lam_R[None, :, None] * dzu)
    return out


def mode_structure_boussinesq(k, U, N, z_eval, hydrostatic=False):
    """w_hat(z) for unit w_hat(0) and lambda (w ~ e^{lambda z})."""
    l2 = (N / U) ** 2
    s = l2 - (0.0 if hydrostatic else k * k)
    if s > 0:
        lam = 1j * np.sign(k * U) * np.sqrt(s)
    else:
        lam = -np.sqrt(-s) + 0j
    return np.exp(lam * np.asarray(z_eval, dtype=float)), lam


# --- full fields -------------------------------------------------------------


def linear_fields(
    x_SI,
    z_SI,
    h_SI,
    U,
    bg,
    regime="compressible",
    hydrostatic=False,
    z_R=24e3,
    dz=5.0,
):
    """Steady linear (u', w, theta', pi') at points (x_i, z_ij).

    x_SI: (Nx,) equispaced periodic sample points (the terrain spectrum is
    the DFT of h_SI on them); z_SI: (Nx, Nz) physical heights. Returns a
    dict of (Nx, Nz) arrays (pi' is the dimensionless Exner perturbation).
    ``bg`` is a :class:`Background`; the Boussinesq regime uses only its
    surface values N(0), theta(0) (and rho(0) for the drag).
    """
    if regime not in REGIMES:
        raise ValueError(f"regime must be one of {REGIMES}")
    x_SI = np.asarray(x_SI, dtype=float)
    z_SI = np.asarray(z_SI, dtype=float)
    nx = x_SI.size
    dx = x_SI[1] - x_SI[0]
    k_all = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    h_hat = np.fft.fft(h_SI) / nx  # h = sum_k h_hat e^{i k (x - x0)}

    # real field: k > 0 modes doubled (the -k mode is the conjugate); the
    # Nyquist mode (even nx) is its own conjugate partner and counted once
    idx = np.arange(1, nx // 2 + 1)
    ks = np.abs(k_all[idx])
    fac = np.where(2 * idx == nx, 1.0, 2.0)
    keep = np.abs(h_hat[idx]) > 1e-14 * np.abs(h_hat).max()
    idx, ks, fac = idx[keep], ks[keep], fac[keep]
    a_w = 1j * ks * U * h_hat[idx]  # w_hat(0) = i k U h_hat

    zu, inv = np.unique(z_SI, return_inverse=True)
    inv = inv.reshape(nx, -1)

    if regime == "boussinesq":
        N0 = np.sqrt(bg.N2(0.0))
        th0 = bg.theta(0.0)
        s = (N0 / U) ** 2 - (0.0 if hydrostatic else 1.0) * ks**2
        lam = np.where(
            s > 0,
            1j * np.sign(ks * U) * np.sqrt(np.abs(s)),
            -np.sqrt(np.abs(s)) + 0j,
        )
        w_hat = a_w[:, None] * np.exp(lam[:, None] * zu[None, :])
        u_hat = (1j * lam / ks)[:, None] * w_hat  # i k u + w_z = 0
        pi_hat = -U * u_hat / (bg.cp * th0)
        theta_hat = 1j * (th0 * N0**2 / bg.g) * w_hat / (ks * U)[:, None]
    else:
        Y = mode_structures_compressible(bg, ks, U, zu, hydrostatic, z_R, dz)
        w_hat = a_w[:, None] * Y[0]
        pi_hat = a_w[:, None] * Y[1]
        u_hat = -bg.cp * bg.theta(zu)[None, :] * pi_hat / U
        theta_hat = 1j * bg.dtheta_dz(zu)[None, :] * w_hat / (ks * U)[:, None]

    # phases from the first sample: e^{i k (x - x0)}  (nk, nx)
    phase = fac[:, None] * np.exp(1j * ks[:, None] * (x_SI - x_SI[0])[None, :])
    out = {}
    for q, qh in (("w", w_hat), ("u", u_hat), ("theta", theta_hat), ("pi", pi_hat)):
        res = np.empty(inv.shape)
        for i in range(nx):  # per column (heights are column-dependent)
            res[i] = np.real(phase[:, i] @ qh[:, inv[i]])
        out[q] = res.reshape(z_SI.shape)
    return out


def drag_profile(x_SI, z_levels, h_SI, U, bg, **kw):
    """D(z) = -int rho_b u' w' dx on flat levels z [N/m]."""
    z_levels = np.asarray(z_levels, dtype=float)
    z = np.broadcast_to(z_levels, (x_SI.size, z_levels.size))
    f = linear_fields(x_SI, z, h_SI, U, bg, **kw)
    dx = x_SI[1] - x_SI[0]
    if kw.get("regime", "compressible") == "boussinesq":
        rho = bg.rho(0.0) + 0.0 * z_levels
    else:
        rho = bg.rho(z_levels)
    return -rho * np.sum(f["u"] * f["w"], axis=0) * dx


def smith_drag(U, N, h0, rho0):
    """Smith (1980) Agnesi drag (pi/4) rho0 N U h0^2 [N/m]."""
    return 0.25 * np.pi * rho0 * N * U * h0**2


# --- model extraction + comparison --------------------------------------------


def _inner(arr, ndim):
    if ndim == 2:
        return arr[2:-2, 2:-2]
    return arr[2:-2, 2:-2, 0]


def sim_fields(mem, ud):
    """Raw model perturbation fields in SI (no anelastic rescaling).

    Cells: u' = u - U, w, theta' = (Y - Y0) T_ref at physical cell heights.
    Nodes: pi' = Msq * p2_nodes (p2 is the perturbation Exner / Msq).
    """
    ndim = mem.elem.ndim
    rho = _inner(mem.sol.rho, ndim)
    u = _inner(mem.sol.rhou, ndim) / rho * ud.u_ref
    w = _inner(mem.sol.rhov, ndim) / rho * ud.u_ref
    Y = _inner(mem.sol.rhoY, ndim) / rho
    Y0 = _inner(mem.npf.HydroState.Y0, ndim)
    x = mem.elem.x[2:-2] * ud.h_ref
    z = _inner(mem.elem.metric.z, ndim) * ud.h_ref
    zeta = mem.elem.y[2:-2] * ud.h_ref
    p2 = mem.npf.p2_nodes
    zn_full = mem.node.metric.z
    if ndim == 3:
        p2 = p2[..., 0]
        zn_full = zn_full[..., 0]
    xn = mem.node.x[2:-2] * ud.h_ref
    return {
        "x": x,
        "z": z,
        "zeta": zeta,
        "u": u - ud.U0,
        "w": w,
        "theta": (Y - Y0) * ud.T_ref,
        "xn": xn,
        "zn": zn_full[2:-2, 2:-2] * ud.h_ref,
        "pi": p2[2:-2, 2:-2] * ud.Msq,
        "rho_SI": rho * ud.p_ref / (ud.R_gas * ud.T_ref),
    }


def oracle_for_sim(sf, ud, regime="compressible", hydrostatic=False, z_R=24e3):
    """Oracle fields at the simulation's cell (and node) points."""
    bg = Background.from_ud(ud)
    h = ud.orography(sf["x"] / ud.h_ref, 0.0) * ud.h_ref
    ref = linear_fields(
        sf["x"], sf["z"], h, ud.U0, bg, regime=regime, hydrostatic=hydrostatic, z_R=z_R
    )
    # nodes: the periodic node row duplicates the first node; drop it for the
    # DFT and evaluate the terrain spectrum on the node points
    xn = sf["xn"][:-1]
    hn = ud.orography(xn / ud.h_ref, 0.0) * ud.h_ref
    refn = linear_fields(
        xn,
        sf["zn"][:-1],
        hn,
        ud.U0,
        bg,
        regime=regime,
        hydrostatic=hydrostatic,
        z_R=z_R,
    )
    ref["pi_nodes"] = refn["pi"]
    return ref


def demean(q):
    return q - q.mean(axis=0, keepdims=True)


def rel_l2(sim, ref, mask):
    s, r = demean(sim), demean(ref)
    return float(
        np.linalg.norm(np.where(mask, s - r, 0.0))
        / np.linalg.norm(np.where(mask, r, 0.0))
    )


def compare(sf, ref, windows, dx_SI, rho_b=None):
    """Rel-L2 per field and window (windows: name -> (z_lo, z_hi) on eta).

    Windows are defined on the eta (computational) level so a window is a
    set of full grid rows; 'z' heights still vary with x near the hill.
    """
    out = {}
    zeta = sf["zeta"]
    for name, (lo, hi) in windows.items():
        rows = (zeta >= lo) & (zeta <= hi)
        mask = np.broadcast_to(rows, sf["w"].shape)
        res = {q: rel_l2(sf[q], ref[q], mask) for q in ("w", "u", "theta")}
        if "pi_nodes" in ref:
            # node rows: mask by node eta, last periodic column dropped
            zn_eta = sf["zn"][-1, :]
            nmask = np.broadcast_to(
                (zn_eta >= lo) & (zn_eta <= hi), ref["pi_nodes"].shape
            )
            res["pi"] = rel_l2(sf["pi"][:-1], ref["pi_nodes"], nmask)
        out[name] = res
    # momentum flux on eta rows (rho-weighted, demeaned u')
    flux_sim = -np.sum(sf["rho_SI"] * demean(sf["u"]) * sf["w"], axis=0) * dx_SI
    flux_ref = -np.sum(sf["rho_SI"] * ref["u"] * ref["w"], axis=0) * dx_SI
    out["flux_sim"], out["flux_ref"] = flux_sim, flux_ref
    return out
