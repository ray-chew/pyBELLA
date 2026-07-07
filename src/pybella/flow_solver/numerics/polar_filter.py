"""FFT-in-longitude polar filter (Stage F, F3).

Near a lat-lon pole the zonal grid spacing a*cos(phi)*dlambda shrinks to
zero, so the high zonal wavenumbers have an effective Courant number > 1
and would force a vanishing dt. The classic fix (Takacs and others, the
FV dynamical-core family) damps exactly those modes with a Fourier filter
in longitude, restoring a sane dt while leaving the well-resolved tropical
flow untouched.

Transfer function (per zonal wavenumber k, latitude phi):

    r(k, phi) = min( 1, [ (cos phi / cos phi_c) / sin(pi k / N) ]^p ),  k >= 1
    r(0, .)   = 1                                (the zonal mean is kept)

with N the interior longitude count, phi_c the onset latitude (default
60 deg) and p the sharpness (default 2). For |phi| <= phi_c the bracket is
>= 1 so r = 1 automatically (no damping in the tropics); poleward it
progressively removes the CFL-violating modes, and at the pole ring only
k = 0 (the ring mean) survives.

Conservation: the filter acts on the J-WEIGHTED conservative fields and
leaves k = 0 untouched, so every longitude ring's J-weighted integral —
hence global mass, momentum, rhoY, tracer — is conserved to roundoff, even
over terrain (J longitude-dependent). Activated by ``ud.polar_filter``
(a :class:`PolarFilter`); ``None`` -> the whole solver is bit-identical.
"""

import numpy as np

from ...utils import axes

_FIELDS = ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX")


class PolarFilter:
    """Polar-filter configuration: onset latitude ``phi_c`` (radians) and
    sharpness exponent ``p``."""

    def __init__(self, phi_c, p=2.0):
        self.phi_c = float(phi_c)
        self.p = float(p)


def transfer(N, cosphi, phi_c, p):
    """Transfer factors ``r[nmodes, nphi]`` for the rfft modes k = 0..N/2.

    ``cosphi`` is cos(latitude) at the interior phi rows; the k = 0 row is
    forced to 1 so the zonal mean (ring integral) is preserved exactly.
    """
    k = np.arange(N // 2 + 1)
    sin_k = np.sin(np.pi * k / N)
    sin_k[0] = 1.0  # placeholder; the k=0 row is overwritten to 1 below
    ratio = (cosphi[None, :] / np.cos(phi_c)) / sin_k[:, None]
    r = np.minimum(1.0, np.abs(ratio) ** p)
    r[0, :] = 1.0
    return r


def cfl_cap(elem, ud):
    """Per-cell longitude signal-speed cap min(1, cos phi / cos phi_c), or
    ``None`` when the filter is inactive.

    With the filter on, the surviving zonal modes have an effective speed
    scaled by ~cos phi / cos phi_c near the pole, so the longitude CFL is
    relieved by exactly this factor — this is what buys the larger dt.
    """
    cfg = getattr(ud, "polar_filter", None)
    m = elem.metric
    if cfg is None or m is None or m.vertical_line:
        return None
    phi_axis = m.cart_haxes[1]
    shape = [1] * elem.ndim
    shape[phi_axis] = -1
    phi = axes.coords_along(elem, phi_axis).reshape(shape)
    return np.minimum(1.0, np.cos(phi) / np.cos(cfg.phi_c))


def apply(mem, ud):
    """Filter the conservative fields in place (canonical orientation)."""
    cfg = getattr(ud, "polar_filter", None)
    if cfg is None:
        return
    elem, sol = mem.elem, mem.sol
    m = elem.metric
    if m is None or m.vertical_line:
        return

    lam_axis = m.cart_haxes[0]
    phi_axis = m.cart_haxes[1]
    v_axis = m.cart_v
    assert lam_axis < phi_axis, "polar filter runs in canonical orientation"

    sc = [int(s) for s in elem.sc]
    igl, igp, igr = (int(elem.igs[a]) for a in (lam_axis, phi_axis, v_axis))
    N = sc[lam_axis] - 2 * igl

    sl = [slice(None)] * elem.ndim
    sl[lam_axis] = slice(igl, sc[lam_axis] - igl)
    sl[phi_axis] = slice(igp, sc[phi_axis] - igp)
    sl[v_axis] = slice(igr, sc[v_axis] - igr)
    sl = tuple(sl)

    phi = axes.coords_along(elem, phi_axis)[igp : sc[phi_axis] - igp]
    r = transfer(N, np.cos(phi), cfg.phi_c, cfg.p)  # (nmodes, nphi)
    # broadcast r over the interior block: modes on lam_axis, r on phi_axis
    bidx = [None] * elem.ndim
    bidx[lam_axis] = slice(None)
    bidx[phi_axis] = slice(None)
    r_bcast = r[tuple(bidx)]

    J = m.J[sl]
    for name in _FIELDS:
        f = getattr(sol, name)
        Jf = J * f[sl]
        F = np.fft.rfft(Jf, axis=lam_axis)
        F *= r_bcast
        f[sl] = np.fft.irfft(F, n=N, axis=lam_axis) / J
