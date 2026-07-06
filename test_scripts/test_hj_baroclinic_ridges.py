"""Hughes & Jablonowski (2023) pt 2 — ridge-triggered baroclinic-wave smoke.

pt 1 (``test_hj_baroclinic``) showed the Ullrich base state on the SMOOTH
shell stays (nearly) steady: the flat background generates only a small
truncation-level meridional adjustment that plateaus at ~0.18 m/s and is NOT
longitude-localised (it sits wherever the discrete residual happens to
concentrate). pt 2 embeds the two H&J midlatitude ridges (Eq. 1) as terrain
geometry (``SphericalTerrainMap``) and re-derives the balanced state on the
tilted coordinate surfaces (the analytic Ullrich state sampled at the
terrain-following height -> the H&J adjusted surface pressure, Eq. 2). That
adjusted state is well-balanced but not PERFECTLY so, and the residual near
the ridges is the intended baroclinic-wave trigger.

This coarse numpy smoke asserts the qualitative dynamical-core behaviour the
paper describes for the FIRST stage of the run: the wave INITIATES FROM THE
RIDGES and nothing blows up. Concretely, over a short (~40 min) coarse run on
the 32 x 12 x 32 shell (the latitude count pinned near 32 so the free-slip
ghost rows past +-80 deg keep J = r^2 cos phi > 0, as in pt 1):

* the state stays finite and the jet stays bounded (no blow-up);
* the analytic momenta are purely zonal at init (no meridional wind);
* the ridges drive a meridional-wind response that is (a) MUCH larger than
  the flat-background pt-1 adjustment (measured ~1.45 m/s, ~8x the 0.18 m/s
  flat plateau) and (b) LONGITUDE-LOCALISED at the ridge centres 72 deg E /
  140 deg E (measured peaks at 73 / 141 deg E). Localisation is the decisive
  signature that the response is ridge-TRIGGERED, not a spurious global
  imbalance -- a broken balance would drive an un-localised O(jet) meridional
  wind, as it does for TC2.

The multi-day maturation into the full Rossby wave train (paper Sect. 4) is
the production device run, pt 3.
"""

import logging

import numpy as np

logging.disable(logging.INFO)

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.flow_solver.utils.boundary import cell_boundary as bdry_c
from pybella.tests import test_hj_baroclinic_ridges as hjr
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.io.debug import NullDebugWriter

#: H&J ridge centre longitudes [deg E] (Table 1)
_RIDGE_LONS = np.array([72.0, 140.0])


def _build(nx, ny, nz, nsteps, dt=None):
    udo = hjr.UserData()
    udo.inx, udo.iny, udo.inz = nx + 1, ny + 1, nz + 1
    udo.stepmax = nsteps
    if dt is not None:
        udo.dtfixed = udo.dtfixed0 = dt
    udo.diag = False
    udo.output_timesteps = False
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = hjr.sol_init(sol, npf, elem, node, th, ud)
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


def test_ridges_initiate_the_wave():
    """The two ridges trigger a growing, ridge-localised meridional wind
    (the baroclinic wave initiating) without blowing up."""
    mem, ud = _build(32, 12, 32, nsteps=8)
    u0, v0 = _winds_ms(mem, ud)
    jet = np.abs(u0).max()
    # the analytic momenta are purely zonal at init, ridges or not
    assert np.abs(v0).max() < 1e-10, np.abs(v0).max()
    assert 25.0 < jet < 30.0, jet  # ~28 m/s Ullrich jet (Fig. 1a)

    mem = _run(mem, ud)
    u1, v1 = _winds_ms(mem, ud)

    for a in ("rho", "rhou", "rhov", "rhow", "rhoY"):
        assert np.all(np.isfinite(getattr(mem.sol, a))), a

    merid_wind = np.abs(v1).max()
    # jet stays bounded: local acceleration over the ridges is modest (~30.5
    # measured, 1.11x); a broken balance would run away
    assert u1.max() < 1.2 * jet, u1.max()

    # (1) the ridges drive meridional wind FAR above the flat-background pt-1
    # adjustment (~0.18 m/s). Measured ~1.45 m/s; gate at 0.5 (~2.8x flat,
    # ~3x under the measured signal) and below O(jet) (not a blow-up)
    assert 0.5 < merid_wind < 0.25 * jet, merid_wind

    # (2) that response is LONGITUDE-LOCALISED at the ridges -- the decisive
    # "ridge-triggered" signature. The longitude of the peak meridional wind
    # (max over r, phi) must sit within 20 deg of a ridge centre (72/140 E);
    # measured 73 / 141 E
    lam_deg = np.degrees(mem.elem.x[2:-2]) % 360.0
    vlon = np.abs(v1).max(axis=(1, 2))
    peak_lon = lam_deg[np.argmax(vlon)]
    dist = np.abs(((peak_lon - _RIDGE_LONS + 180.0) % 360.0) - 180.0).min()
    assert dist < 20.0, (peak_lon, dist)
