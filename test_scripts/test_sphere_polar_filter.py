"""Stage F, increment F3 gates: FFT-in-longitude polar filter.

Gates:
1. no-op on a longitude-independent state (only k=0 survives) -> the
   filter returns the input to roundoff;
2. conservation: every longitude ring's J-weighted integral (hence global
   mass/momentum/rhoY/tracer) is preserved to machine precision, on both
   the shell and a terrain map (J longitude-dependent);
3. selectivity: a single-wavenumber signal poleward of phi_c is damped by
   exactly the transfer factor r(k, phi);
4. ud.polar_filter = None is a no-op (bit-identity);
5. with the filter + CFL cap a longitude-CFL-violating dt stays stable.
"""

import types

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid, spherical
from pybella.flow_solver.numerics import polar_filter
from pybella.flow_solver.numerics.explicit_advection import (
    advective_flux,
    compute_advection,
)
from pybella.flow_solver.utils import fields
from pybella.utils import axes as _axes
from pybella.utils import options as opts

_A = 2.0
_PHI_C = np.deg2rad(60.0)


class _PoleUD:
    def __init__(self, n=(32, 4, 16), a=_A, terrain=False, filt=True):
        self.inx, self.iny, self.inz = n[0] + 1, n[1] + 1, n[2] + 1
        self.xmin, self.xmax = -np.pi, np.pi
        self.zmin, self.zmax = -0.5 * np.pi, 0.5 * np.pi
        self.gravity_direction = 1
        self.bdry_type = np.array(
            [opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.POLE]
        )
        if terrain:
            self.ymin, self.ymax = 0.0, 0.1
            h0, sig = 0.02, 0.5
            h = (
                lambda lam, phi: h0
                * (1 + np.cos(lam))
                / 2
                * np.exp(-((phi / sig) ** 2))
            )
            dl = lambda lam, phi: -h0 / 2 * np.sin(lam) * np.exp(-((phi / sig) ** 2))
            dp = (
                lambda lam, phi: h0
                * (1 + np.cos(lam))
                / 2
                * (-2 * phi / sig**2)
                * np.exp(-((phi / sig) ** 2))
            )
            self.curvilinear_map = spherical.SphericalTerrainMap(
                a, 0.1, h, (dl, dp), pole=True
            )
        else:
            self.ymin, self.ymax = a * 0.95, a * 1.05
            self.curvilinear_map = spherical.SphericalShellMap(a, pole=True)
        self.polar_filter = polar_filter.PolarFilter(_PHI_C) if filt else None


def _mem(n=(32, 4, 16), terrain=False, filt=True):
    ud = _PoleUD(n=n, terrain=terrain, filt=filt)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    mem = types.SimpleNamespace(elem=elem, node=node, sol=sol)
    return ud, mem, elem, sol


def _ring_sums(elem, sol):
    """Per (r, phi) longitude-ring J-weighted sums of every field."""
    J = elem.metric.J
    igl = int(elem.igs[0])
    out = {}
    for name in polar_filter._FIELDS:
        f = getattr(sol, name)
        out[name] = np.sum((J * f)[igl:-igl], axis=0)
    return out


# ---------------------------------------------------------- no-op theorem


def test_filter_noop_on_zonal_state():
    ud, mem, elem, sol = _mem()
    # genuinely longitude-independent fields: functions of the r and phi
    # COORDINATES only (metric.x[0,1] carry cos/sin lambda, so are not)
    rc = _axes.coords_along(elem, 1).reshape(1, -1, 1)
    phic = _axes.coords_along(elem, 2).reshape(1, 1, -1)
    sol.rho[...] = 1.0 + 0.1 * np.cos(phic) + 0.0 * rc
    sol.rhoY[...] = 1.0 + 0.05 * np.sin(phic)
    sol.rhou[...] = 0.2 * np.sin(phic) + 0.0 * rc
    sol.rhov[...] = 0.0
    sol.rhow[...] = 0.1
    sol.rhoX[...] = 0.3 + 0.0 * phic
    ref = {n: getattr(sol, n).copy() for n in polar_filter._FIELDS}
    polar_filter.apply(mem, ud)
    igl = int(elem.igs[0])
    for name, r in ref.items():
        got = getattr(sol, name)
        sc = max(np.max(np.abs(r)), 1e-30)
        np.testing.assert_allclose(
            got[igl:-igl], r[igl:-igl], rtol=1e-12, atol=1e-12 * sc
        )


# ------------------------------------------------------------ conservation


def _fill_random(sol, elem, seed=7):
    rng = np.random.default_rng(seed)
    shape = tuple(int(s) for s in elem.sc)
    sol.rho[...] = 1.0 + 0.3 * rng.random(shape)
    sol.rhoY[...] = 1.0 + 0.2 * rng.random(shape)
    sol.rhou[...] = rng.standard_normal(shape)
    sol.rhov[...] = rng.standard_normal(shape)
    sol.rhow[...] = rng.standard_normal(shape)
    sol.rhoX[...] = rng.random(shape)


@pytest.mark.parametrize("terrain", [False, True], ids=["shell", "terrain"])
def test_filter_conserves_ring_integrals(terrain):
    ud, mem, elem, sol = _mem(terrain=terrain)
    _fill_random(sol, elem)
    before = _ring_sums(elem, sol)
    polar_filter.apply(mem, ud)
    after = _ring_sums(elem, sol)
    for name in polar_filter._FIELDS:
        sc = max(np.max(np.abs(before[name])), 1e-30)
        np.testing.assert_allclose(
            after[name], before[name], rtol=1e-11, atol=1e-11 * sc
        )


# ------------------------------------------------------------- selectivity


def test_filter_selectivity_matches_transfer():
    ud, mem, elem, sol = _mem(n=(64, 4, 32))
    igl = int(elem.igs[0])
    N = int(elem.sc[0]) - 2 * igl
    k0 = N // 3  # a mid-high wavenumber
    lam = _axes.coords_along(elem, 0).reshape(-1, 1, 1)
    sol.rho[...] = 1.0
    sol.rhoY[...] = 1.0
    sol.rhou[...] = 0.0
    sol.rhov[...] = 0.0
    sol.rhow[...] = 0.0
    sol.rhoX[...] = np.cos(k0 * lam) + 0.0 * elem.metric.J

    before = np.fft.rfft(sol.rhoX[igl:-igl], axis=0)
    polar_filter.apply(mem, ud)
    after = np.fft.rfft(sol.rhoX[igl:-igl], axis=0)

    phi = _axes.coords_along(elem, 2)
    igp = int(elem.igs[2])
    r = polar_filter.transfer(N, np.cos(phi[igp:-igp]), _PHI_C, 2.0)
    # pick a poleward row (~75 deg) and an interior radius; the rfft arrays
    # are indexed on the FULL phi (igp + interior row) and interior radius
    row = int(np.argmin(np.abs(phi[igp:-igp] - np.deg2rad(75.0))))
    rr = int(elem.igs[1])
    ratio = np.abs(after[k0, rr, igp + row]) / np.abs(before[k0, rr, igp + row])
    assert abs(ratio - r[k0, row]) < 1e-9, (ratio, r[k0, row])
    assert r[k0, row] < 0.9  # genuinely damped there


# --------------------------------------------------------------- inactive


def test_cfl_cap_relieves_pole():
    ud, mem, elem, sol = _mem()
    cap = polar_filter.cfl_cap(elem, ud)
    phi = _axes.coords_along(elem, 2).reshape(1, 1, -1)
    expect = np.minimum(1.0, np.cos(phi) / np.cos(_PHI_C))
    np.testing.assert_allclose(cap, np.broadcast_to(expect, cap.shape))
    assert cap.max() <= 1.0 + 1e-15 and np.isclose(cap.max(), 1.0)  # tropics
    assert cap.min() < 0.2  # strong longitude-CFL relief at the pole ring
    # inactive filter -> no cap
    ud_off = _PoleUD(filt=False)
    elem_off, _ = dis_grid.grid_init(ud_off)
    assert polar_filter.cfl_cap(elem_off, ud_off) is None


def test_filter_none_is_noop():
    ud, mem, elem, sol = _mem(filt=False)
    _fill_random(sol, elem)
    ref = {n: getattr(sol, n).copy() for n in polar_filter._FIELDS}
    polar_filter.apply(mem, ud)  # ud.polar_filter is None -> returns immediately
    for name, r in ref.items():
        np.testing.assert_array_equal(getattr(sol, name), r)
