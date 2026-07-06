"""Stage-D5 physics gates: nonhydrostatic gravity wave on the small-planet
compressible shell (``pybella.tests.test_sphere_gw``).

The gates prove the general (curved-metric) radial-gravity chain is both
CONSERVATIVE and second-order-accurate, and that its wave physics matches
the planar Baldauf & Brdar (2013) linear oracle in the large-radius limit:

- ``test_gravity_wave_conservation`` — J-weighted mass and P = rho*Y are
  conserved to machine precision under dynamics. This exercises the
  well-balanced general free-slip wall: reflecting only the up-velocity
  leaves an O(dz^2) mass leak from the hydrostatic rho/rhoY jump across the
  radial walls; making the rhoY flux Y*(N_v.m) exactly odd (image-cell
  N_v.e_up) cancels the wall flux to roundoff (see
  ``cell_boundary._calculate_ghost_values``).
- ``test_gravity_wave_self_convergence`` — the radial-velocity field
  converges at ~2nd order under joint (lambda, r) refinement (Richardson
  on the fields interpolated to the coarse grid).
- ``test_large_radius_limit_baldauf_brdar`` — at increasing planet radius
  the near-equator lambda-height slice converges to the planar B&B linear
  solution of its own initial condition (reusing
  ``baldauf_brdar_analytic.evolve_linear``); the mismatch shrinks with a.
"""

import logging

import numpy as np
import pytest

logging.disable(logging.INFO)

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import spherical, time_update
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.tests import test_sphere_gw as gw
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState
from pybella.utils.io.debug import NullDebugWriter


def _build(nx, ny, nz, nsteps, X=125.0, dt=None, amp=None):
    udo = gw.UserData()
    udo.X = X
    udo.planet_radius = gw._A_EARTH / X
    a_nd = udo.planet_radius / udo.h_ref
    depth = udo.depth_m / udo.h_ref
    udo.ymin, udo.ymax = a_nd, a_nd + depth
    udo.curvilinear_map = spherical.SphericalShellMap(a_nd, frozen_radius=False)
    udo.inx, udo.iny, udo.inz = nx + 1, ny + 1, nz + 1
    udo.stepmax = nsteps
    if dt is not None:
        udo.dtfixed = udo.dtfixed0 = dt
    if amp is not None:
        udo.pert_amplitude = amp
    udo.diag = False
    udo.output_timesteps = False

    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = gw.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    return mem, ud


def _run(mem, ud):
    mem = time_update.do(
        mem, ud, tout=1.0e9, bld=None, writer=None, debug_writer=NullDebugWriter()
    )
    assert mem.time.step == ud.stepmax
    return mem


# --------------------------------------------------------------- conservation


@pytest.mark.parametrize("ny", [8, 16], ids=["ny8", "ny16"])
def test_gravity_wave_conservation(ny):
    """J-weighted mass and P are conserved to roundoff by the well-balanced
    general free-slip wall, AND a gravity wave of sane magnitude develops."""
    mem, ud = _build(64, ny, 8, nsteps=40)
    inner = (slice(2, -2),) * 3
    J = mem.elem.metric.J[inner]
    mass0 = np.sum(J * mem.sol.rho[inner])
    P0 = np.sum(J * mem.sol.rhoY[inner])

    mem = _run(mem, ud)

    mass1 = np.sum(J * mem.sol.rho[inner])
    P1 = np.sum(J * mem.sol.rhoY[inner])
    assert abs(mass1 - mass0) / mass0 < 1e-10, (mass1 - mass0) / mass0
    assert abs(P1 - P0) / P0 < 1e-10, (P1 - P0) / P0

    # a nonhydrostatic gravity-wave response actually developed (not at rest,
    # not blown up): radial velocity in a physically sane band
    w = np.max(np.abs(mem.sol.rhov[inner] / mem.sol.rho[inner])) * ud.u_ref
    assert 1e-4 < w < 5.0, f"gravity-wave w response out of range: {w:.3e} m/s"


# ------------------------------------------------------------ self-convergence


def _run_field(nx, T):
    """Radial-momentum field + coords after physical time T, refining
    (lambda, r) together with dt ~ 1/nx (fixed-time joint refinement)."""
    ny = max(4, nx // 8)
    dt = 0.006833 * 64.0 / nx
    nsteps = int(round(T / dt))
    mem, ud = _build(nx, ny, 8, nsteps=nsteps, dt=dt)
    mem = _run(mem, ud)
    inner = (slice(2, -2),) * 3
    w = (mem.sol.rhov / mem.sol.rho)[inner]
    return mem.elem.x[2:-2], mem.elem.y[2:-2], mem.elem.z[2:-2], w


def test_gravity_wave_self_convergence():
    """Second-order self-convergence of the developed wave under joint
    (lambda, r) refinement: ||w_n - w_2n|| drops ~4x per doubling
    (Richardson on the fields interpolated to the coarse grid). The curved
    metric introduces no order reduction."""
    from scipy.interpolate import RegularGridInterpolator as RGI

    # (32, 64, 128) / T = 1.5: the n=32 bump (~2 cells) is marginal but the
    # 64-128 pair is well resolved -> clean asymptotic ratio (calibrated 3.92).
    # Coarser/shorter (24/48/96, T=1.2) under-resolves the bump (ratio ~2.2).
    T = 1.5
    ns = (32, 64, 128)
    res = {n: _run_field(n, T) for n in ns}

    def interp(src, tgt):
        ls, rs, ps, _ = res[src]
        lt, rt, pt, _ = res[tgt]
        f = RGI((ls, rs, ps), res[src][3], bounds_error=False, fill_value=None)
        L, R, P = np.meshgrid(lt, rt, pt, indexing="ij")
        return f((L, R, P))

    d1 = np.sqrt(np.mean((res[ns[0]][3] - interp(ns[1], ns[0])) ** 2))
    d2 = np.sqrt(np.mean((res[ns[1]][3] - interp(ns[2], ns[1])) ** 2))
    ratio = d1 / d2
    # 4 = exact 2nd order; gate well above 1st order (2), calibrated ~3.9
    assert ratio > 3.0, f"self-convergence ratio {ratio:.2f} (expect ~4)"


# ----------------------------------------------------- B&B large-radius limit

# The gravity-wave case shares the isothermal Baldauf & Brdar background
# (T_ref = 250 K, u_ref = 10 m/s, h_ref = R T_ref / g), so the near-equator
# (lambda, r) slice of the sphere run can be compared to the SAME planar
# linear oracle (baldauf_brdar_analytic) used by test_igw_analytic, with
# x = a * lambda. Fixing the PHYSICAL perturbation (Lambda ~ 1/a) and the
# physical resolution (nx ~ a) as the planet radius grows, the slice
# converges to the planar solution (curvature O((L/a)^2) and acoustic
# wrap-around O(1/a) both vanish).
_H_REF = 287.05 * 250.0 / 9.80665
_A0 = gw._A_EARTH / 125.0 / _H_REF  # a_nd at X = 125
_A_PHYS = _A0 * 0.4  # fixed physical perturbation half-width [h_ref]
_PHI_BAND = 0.12  # narrow equatorial band [rad]


def _run_bb(X, T=0.8, nz=12, nx0=64):
    from pybella.tests import baldauf_brdar_analytic as bb

    a_nd = gw._A_EARTH / X / _H_REF
    nx = int(round(nx0 * a_nd / _A0 / 4) * 4)
    dt = 0.006833 * a_nd / _A0 * nx0 / nx  # ~ const across a (dx_phys, CFL)
    udo = gw.UserData()
    udo.X = X
    udo.planet_radius = gw._A_EARTH / X
    depth = udo.depth_m / udo.h_ref
    udo.ymin, udo.ymax = a_nd, a_nd + depth
    udo.curvilinear_map = spherical.SphericalShellMap(a_nd, frozen_radius=False)
    udo.pert_halfwidth = _A_PHYS / a_nd
    udo.phi_band = _PHI_BAND
    udo.zmin, udo.zmax = -_PHI_BAND, _PHI_BAND
    udo.inx, udo.iny, udo.inz = nx + 1, nz + 1, 4 + 1
    udo.dtfixed = udo.dtfixed0 = dt
    udo.stepmax = int(round(T / dt))
    udo.diag = False
    udo.output_timesteps = False
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    sol = gw.sol_init(sol, npf, elem, node, th, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    par = bb.IGWParams(ud)
    par.f = 0.0  # nonrotating

    def slice_SI(m):
        e = m.elem
        phi = e.z[2:-2]
        jz = int(np.argmin(np.abs(phi)))
        phic = phi[jz]
        sl = (slice(2, -2), slice(2, -2), 2 + jz)
        rho = m.sol.rho[sl]
        rhoY = m.sol.rhoY[sl]
        ru, rv, rw = m.sol.rhou[sl], m.sol.rhov[sl], m.sol.rhow[sl]
        lam = e.x[2:-2].reshape(-1, 1)
        cl, sll = np.cos(lam), np.sin(lam)
        cp, sp = np.cos(phic), np.sin(phic)
        el = (-sll, cl, 0.0 * lam)  # e_lambda (zonal)
        er = (cp * cl, cp * sll, -sp + 0 * lam)  # e_r (vertical)
        ep = (-sp * cl, -sp * sll, -cp + 0 * lam)  # e_phi (out of plane)
        pr = lambda ev: (ev[0] * ru + ev[1] * rv + ev[2] * rw) / rho
        z = (e.y[2:-2] - a_nd) * ud.h_ref
        return {
            "u": pr(el) * ud.u_ref,
            "vo": pr(ep) * ud.u_ref,
            "w": pr(er) * ud.u_ref,
            "p": rhoY**par.gamma * par.p_s - par.p0(z)[None, :],
            "rho": rho * par.rho_ref - (par.p0(z) / (par.R * par.T0))[None, :],
        }, z

    ic, z = slice_SI(mem)
    mem = _run(mem, ud)
    end, _ = slice_SI(mem)
    x_len = 2.0 * np.pi * a_nd * ud.h_ref
    ref, diag = bb.evolve_linear(ic, x_len, z, mem.time.t * ud.t_ref, par, refine=2)

    dm = lambda q: q - q.mean(axis=0, keepdims=True)
    errs = {}
    for k in ("u", "w", "p", "rho"):
        s, r = dm(end[k]), ref[k]
        errs[k] = np.linalg.norm(s - r) / np.linalg.norm(r)
    amp = np.abs(dm(end["u"])).max() / np.abs(ref["u"]).max()
    return {"nx": nx, "errs": errs, "amp": amp, "edrift": diag["energy_drift"]}


def test_large_radius_limit_baldauf_brdar():
    """The near-equator slice matches the planar Baldauf & Brdar linear
    oracle and CONVERGES to it as the planet radius grows.

    The zonal wave field u is the clean indicator: at X=125 (a~7 scale
    heights) its rel-L2 to the oracle already sits inside the planar case's
    own gate (test_igw_analytic uses 0.40; pyBELLA-vs-exact-linear carries
    an intrinsic ~0.2-0.3 from semi-implicit acoustics + O((omega dt)^2)
    phase error), with amplitude ratio ~1. Halving the planet (a doubles)
    halves the error (~1/a: measured 0.187 -> 0.098 -> 0.051 for
    X=125,62.5,31.25). p/rho are only loosely bounded (the semi-implicit
    scheme damps the acoustic pressure the exact-linear oracle keeps); w is
    a 10x-smaller field, left ungated like the planar oracle's loose fields.
    """
    r1 = _run_bb(125.0)  # small planet (worst case)
    r2 = _run_bb(62.5)  # 2x radius

    assert r1["edrift"] < 1e-9, r1["edrift"]  # comparator-internal exactness
    # zonal wave field matches the oracle at the planar-case level
    assert r1["errs"]["u"] < 0.40, r1["errs"]["u"]
    assert 0.7 < r1["amp"] < 1.3, r1["amp"]
    # loose sanity on the density wave (planar gate 0.50)
    assert r1["errs"]["rho"] < 0.55, r1["errs"]["rho"]
    # large-radius limit: the zonal-field error shrinks toward the plane
    assert r2["errs"]["u"] < 0.8 * r1["errs"]["u"], (r1["errs"]["u"], r2["errs"]["u"])
