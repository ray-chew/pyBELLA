"""Hughes & Jablonowski (2023) on the FULL pole-to-pole sphere (Stage F, F8).

The global (``BdryType.POLE`` + polar filter) analogues of the +-80 deg
channel gates ``test_hj_baroclinic`` (flat-background steadiness) and
``test_hj_baroclinic_ridges`` (ridge-triggered initiation). Same qualitative
dynamical-core assertions; the only differences are the pole-to-pole domain
(the channel's "phi >= 32" and 89.5 deg clamp are gone — see the case
docstrings) and the filter (phi_c = 70 deg, poleward of the 45 deg N ridges).

* flat background: the balanced Ullrich state barely moves (jet bounded, zonal
  drift tiny, meridional wind stays a small adjustment). The filter is a near
  no-op — the background is longitude-independent (k = 0), which the filter
  keeps exactly.
* ridges: the two midlatitude ridges drive a meridional-wind response FAR
  above the flat adjustment, LONGITUDE-LOCALISED at 72 / 140 deg E — the
  ridge-triggered signature, without blowing up.
* device: a short jax-device-vs-numpy reproduction (F7b) — the ridge-global
  case is the union of already-validated device paths + the poles.

Coarse numpy smoke; the production multi-day device run stays GPU-gated (pt 3
protocol). Skips the device gate cleanly when jax is not installed.
"""

import importlib.util
import logging

import numpy as np
import pytest

logging.disable(logging.INFO)

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import test_hj_baroclinic_global as hjg
from pybella.tests import test_hj_baroclinic_ridges_global as hjrg
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.io.debug import NullDebugWriter

#: H&J ridge centre longitudes [deg E] (Table 1)
_RIDGE_LONS = np.array([72.0, 140.0])


def _build(case, nx, ny, nz, nsteps, backend="numpy"):
    udo = case.UserData()
    udo.inx, udo.iny, udo.inz = nx + 1, ny + 1, nz + 1
    udo.stepmax = nsteps
    udo.diag = False
    udo.output_timesteps = False
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    ud.backend = backend
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = case.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    bdry_c.set_ghost_cells(mem, ud)
    return mem, ud


def _winds_ms(mem, ud):
    """Interior zonal (e_lambda) and meridional (e_phi) wind fields [m/s]."""
    e = mem.elem
    i2 = (slice(2, -2), slice(2, -2), slice(2, -2))
    lam = e.x[2:-2].reshape(-1, 1, 1)
    phi = e.z[2:-2].reshape(1, 1, -1)
    cl, sl = np.cos(lam), np.sin(lam)
    cp, sp = np.cos(phi), np.sin(phi)
    e_lam = (-sl, cl, 0.0 * lam + 0.0 * phi)
    e_phi = (-sp * cl, -sp * sl, -cp + 0.0 * lam)
    rho = mem.sol.rho[i2]
    ru, rv, rw = mem.sol.rhou[i2], mem.sol.rhov[i2], mem.sol.rhow[i2]
    proj = lambda ev: (ev[0] * ru + ev[1] * rv + ev[2] * rw) / rho
    return proj(e_lam) * ud.u_ref, proj(e_phi) * ud.u_ref


def _run(mem, ud):
    mem = time_update.do(
        mem, ud, tout=1.0e9, bld=None, writer=None, debug_writer=NullDebugWriter()
    )
    assert mem.time.step == ud.stepmax
    return mem


def test_global_background_stays_steady():
    """No topography -> the balanced Ullrich background barely moves on the
    full pole-to-pole sphere (the pt-1 tripwire, global)."""
    mem, ud = _build(hjg, 32, 12, 32, nsteps=8)
    u0, v0 = _winds_ms(mem, ud)
    jet = np.abs(u0).max()
    assert np.abs(v0).max() < 1e-10, np.abs(v0).max()  # purely zonal at init
    assert 25.0 < jet < 30.0, jet  # ~28 m/s Ullrich jet (Fig. 1a)

    mem = _run(mem, ud)
    u1, v1 = _winds_ms(mem, ud)
    for a in ("rho", "rhou", "rhov", "rhow", "rhoY"):
        assert np.all(np.isfinite(getattr(mem.sol, a))), a

    zonal_drift = np.sqrt(np.mean((u1 - u0) ** 2))
    merid_wind = np.abs(v1).max()
    assert u1.max() < 1.05 * jet, u1.max()  # jet bounded (no blow-up)
    assert zonal_drift < 0.01 * jet, zonal_drift  # zonal field barely drifts
    # balance tripwire: the adjustment stays a small fraction of the jet (a
    # broken balance / wrong Coriolis would drive O(jet), as for TC2)
    assert merid_wind < 0.05 * jet, merid_wind


def test_global_ridges_initiate_the_wave():
    """The two ridges trigger a growing, ridge-localised meridional wind on the
    full pole-to-pole sphere (the pt-2 signature, global)."""
    mem, ud = _build(hjrg, 32, 12, 32, nsteps=8)
    u0, v0 = _winds_ms(mem, ud)
    jet = np.abs(u0).max()
    assert np.abs(v0).max() < 1e-10, np.abs(v0).max()  # purely zonal at init
    assert 25.0 < jet < 30.0, jet

    mem = _run(mem, ud)
    u1, v1 = _winds_ms(mem, ud)
    for a in ("rho", "rhou", "rhov", "rhow", "rhoY"):
        assert np.all(np.isfinite(getattr(mem.sol, a))), a

    merid_wind = np.abs(v1).max()
    assert u1.max() < 1.2 * jet, u1.max()  # jet bounded
    # (1) meridional wind FAR above the flat-background adjustment, below O(jet)
    assert 0.5 < merid_wind < 0.25 * jet, merid_wind
    # (2) LONGITUDE-LOCALISED at the ridges (the ridge-triggered signature):
    # the peak meridional-wind longitude sits within 20 deg of a ridge centre
    lam_deg = np.degrees(mem.elem.x[2:-2]) % 360.0
    vlon = np.abs(v1).max(axis=(1, 2))
    peak_lon = lam_deg[np.argmax(vlon)]
    dist = np.abs(((peak_lon - _RIDGE_LONS + 180.0) % 360.0) - 180.0).min()
    assert dist < 20.0, (peak_lon, dist)


@pytest.mark.skipif(importlib.util.find_spec("jax") is None, reason="jax not installed")
def test_ridges_global_device_reproduces_numpy():
    """The device-resident JAX backend reproduces numpy for the pole-to-pole
    ridge case (F8 device gate). This is the union of already-validated device
    paths — the general sphere metric with radial gravity + field-mode
    compressible HydroState (channel ridges) and the pole machinery (F7b: pole
    ghost fold, elliptic collapse, polar filter in the step loop). Compressible
    with no projection, so it never hits the field-mode+incompressible guard.

    Agreement is at the per-step bicgstab Krylov floor, which for THIS case
    sits higher than any other sphere gate: the deep compressible shell + the
    terrain tilt + the pole-ring collapse together make the elliptic system the
    worst-conditioned in the suite (dev_notes/sphere_poles_plan.md F8), so the
    two backends' ulp-different systems select Krylov members further apart,
    and XLA's run-to-run CPU reduction-order variation makes the gap fluctuate.
    Measured (3 steps, 32x12x32): scalars ~9e-6, momenta/rho ~1e-4, p2 relative
    ~1.4e-4. Thresholds sit ~3-5x above that floor and ~100x below the ~1e-2 a
    genuinely broken device path produces (e.g. a missing pole ghost branch) —
    the machine-precision reproduction is gated elsewhere (the window and TC2
    device gates), this one certifies the union path stays at the floor and
    finite."""
    n = 3
    mem_np, ud_np = _build(hjrg, 32, 12, 32, nsteps=n)
    mem_dv, ud_dv = _build(hjrg, 32, 12, 32, nsteps=n, backend="jax-device")
    mem_np = _run(mem_np, ud_np)
    mem_dv = _run(mem_dv, ud_dv)
    assert getattr(mem_dv, "_device_compile_count", None) == 2  # parity 0/1

    inner = (slice(2, -2), slice(2, -2), slice(2, -2))
    rho = np.asarray(mem_np.sol.rho)[inner]
    for name in ("rho", "rhoY", "rhoX"):
        d = float(
            np.max(
                np.abs(
                    np.asarray(getattr(mem_dv.sol, name))[inner]
                    - np.asarray(getattr(mem_np.sol, name))[inner]
                )
            )
        )
        assert d < 5e-5, f"{name} (absolute): {d:.3e}"
    for name in ("rhou", "rhov", "rhow"):
        d = float(
            np.max(
                np.abs(
                    (
                        np.asarray(getattr(mem_dv.sol, name))[inner]
                        - np.asarray(getattr(mem_np.sol, name))[inner]
                    )
                    / rho
                )
            )
        )
        assert d < 5e-4, f"{name} (momentum/rho): {d:.3e}"
    p2_np = np.asarray(mem_np.npf.p2_nodes)
    dp2 = float(np.max(np.abs(np.asarray(mem_dv.npf.p2_nodes) - p2_np)))
    assert dp2 < 5e-4 * np.abs(p2_np).max(), f"p2_nodes (relative): {dp2:.3e}"
    for name in ("rho", "rhoY", "rhou", "rhov", "rhow"):
        assert np.all(np.isfinite(np.asarray(getattr(mem_dv.sol, name)))), name
