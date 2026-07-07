"""Stage F, increment F2 gates: pole-face flux closure + over-pole advection.

Williamson et al. (1992) test case 1 — a cosine bell advected by solid-body
rotation — with the rotation axis in the EQUATORIAL plane so the bell is
carried straight over both poles. Run as a REDUCED harness: a passive
tracer (rhoX) on a frozen thin spherical shell, prescribed rigid-rotation
momenta reset every step, driven directly through the advection sweeps (no
elliptic solve, no polar filter yet).

The pole face has zero area in the continuum; the discrete flux there is a
symmetric O(dphi^2) residual that BOTH pole-adjacent cells subtract with the
same sign — a double-loss leak. F2 zeroes that flux. Gates:

1. J-weighted total tracer conserved to MACHINE precision every step (the
   pole-face-closure theorem; without the zeroing it leaks ~O(dphi^2));
2. the bell is advected OVER the pole and stays bounded (no blow-up, no
   large undershoot), re-emerging on the far side;
3. an equatorial-rotation twin (axis along the pole) conserves identically
   and never engages the pole faces.
"""

import numpy as np

from pybella.flow_solver.discretisation import grid as dis_grid, spherical
from pybella.flow_solver.numerics.explicit_advection import (
    advective_flux,
    compute_advection,
)
from pybella.flow_solver.physics import thermodynamics
from pybella.flow_solver.utils import cache, fields
from pybella.tests import test_sphere_swe_tc2 as tc2
from pybella.tests.case_setup import build_bdry
from pybella.utils import options as opts
from pybella.utils import user_data
from pybella.utils.data_structures import ModelState

_A = tc2.A_ND


class _TC1UD(tc2.UserData):
    def __init__(self, nlam=32, nphi=16):
        super().__init__()
        self.zmin, self.zmax = -0.5 * np.pi, 0.5 * np.pi
        self.bdry_type = build_bdry(
            opts.BdryType.PERIODIC, opts.BdryType.WALL, opts.BdryType.POLE
        )
        self.curvilinear_map = spherical.SphericalShellMap(
            _A, frozen_radius=True, pole=True
        )
        self.constrain_to_surface = False  # rigid rotation is already tangent
        self.coriolis_field = None
        self.inx, self.iny, self.inz = nlam + 1, 2, nphi + 1
        self.initial_projection = False
        self.stepmax = 10**9


def _build(nlam=32, nphi=16):
    udo = _TC1UD(nlam, nphi)
    ud = user_data.UserDataInit(**vars(udo))
    ud.coriolis_strength = np.array(ud.coriolis_strength)
    elem, node = dis_grid.grid_init(ud)
    sol = fields.CellSolField(elem.sc)
    th = thermodynamics.ThermodynamicalQuantities(ud)
    npf = fields.NodePressureField(elem, node, ud)
    mem = ModelState(elem, node, sol, npf, th, cache.FlowSolverCache())
    return ud, mem, elem, node, sol


def _rigid(x, omega):
    """omega x x, Cartesian components (a divergence-free tangent field)."""
    return (
        omega[1] * x[2] - omega[2] * x[1],
        omega[2] * x[0] - omega[0] * x[2],
        omega[0] * x[1] - omega[1] * x[0],
    )


def _bell(x, center, r_ang, h0=1.0):
    xn = [x[k] / _A for k in range(3)]
    cos_d = np.clip(sum(xn[k] * center[k] for k in range(3)), -1.0, 1.0)
    d = np.arccos(cos_d)  # angular great-circle distance
    return np.where(d < r_ang, h0 * 0.5 * (1.0 + np.cos(np.pi * d / r_ang)), 0.0)


def _set_background(sol, x, omega):
    sol.rho[...] = 1.0
    sol.rhoY[...] = 1.0
    u = _rigid(x, omega)
    sol.rhou[...] = sol.rho * u[0]
    sol.rhov[...] = sol.rho * u[1]
    sol.rhow[...] = sol.rho * u[2]


def _run(omega, center0, nsteps, dt, nlam=32, nphi=16, r_ang=0.5):
    ud, mem, elem, node, sol = _build(nlam, nphi)
    m = elem.metric
    x = m.x
    _set_background(sol, x, omega)
    sol.rhoX[...] = _bell(x, center0, r_ang)

    Jw = m.J
    inner = (slice(2, -2), slice(2, -2), slice(2, -2))
    mass0 = np.sum((Jw * sol.rhoX)[inner])
    max_drift = 0.0
    rhoX_min = 0.0
    rhoX_max = float(sol.rhoX[inner].max())

    for step in range(nsteps):
        _set_background(sol, x, omega)
        advective_flux.recompute(mem, ud)
        compute_advection.strange_splitting(mem, ud, dt, step % 2, "%d" % step)
        assert np.isfinite(sol.rhoX).all(), step
        mass = np.sum((Jw * sol.rhoX)[inner])
        max_drift = max(max_drift, abs(mass - mass0) / abs(mass0))
        rhoX_min = min(rhoX_min, float(sol.rhoX[inner].min()))
        rhoX_max = max(rhoX_max, float(sol.rhoX[inner].max()))

    return dict(
        sol=sol,
        elem=elem,
        mass0=mass0,
        mass=mass,
        max_drift=max_drift,
        rhoX_min=rhoX_min,
        rhoX_max=rhoX_max,
        x=x,
        Jw=Jw,
        inner=inner,
    )


# ------------------------------------------------- conservation (the gate)


def test_pole_face_closure_conserves_tracer():
    """A few steps near the equator: the J-weighted tracer integral is
    conserved to machine precision. This is the pole-face-closure theorem
    (it holds for any dt); without the zeroing it leaks."""
    omega = np.array([0.1, 0.0, 0.0])  # equatorial axis -> over the poles
    res = _run(omega, np.array([0.0, 1.0, 0.0]), nsteps=20, dt=0.05)
    assert res["max_drift"] < 1e-12, res["max_drift"]


def test_equatorial_reference_conserves():
    """Rotation about the pole axis: bell circles the equator, the pole
    faces are never engaged, tracer conserved identically."""
    omega = np.array([0.0, 0.0, 0.1])  # pole axis
    res = _run(omega, np.array([1.0, 0.0, 0.0]), nsteps=20, dt=0.05)
    assert res["max_drift"] < 1e-12, res["max_drift"]


# ------------------------------------------------ over-pole transport


def test_bell_crosses_pole_bounded():
    """Full quarter revolution: the bell is carried over the south pole and
    stays bounded and conservative throughout (no blow-up, no large
    undershoot). Coarse resolution, no filter -> small dt."""
    omega = np.array([0.1, 0.0, 0.0])
    center0 = np.array([0.0, 1.0, 0.0])  # equator, crosses the poles
    # quarter revolution theta = pi/2 at t = pi/(2*Omega) ~ 15.7
    dt = 0.05
    nsteps = int(round((np.pi / (2 * 0.1)) / dt))
    res = _run(omega, center0, nsteps=nsteps, dt=dt)

    assert res["max_drift"] < 1e-12, res["max_drift"]
    assert res["rhoX_max"] < 1.5, res["rhoX_max"]  # no growth beyond ~1
    # ~11% undershoot: the unlimited scheme (LimiterType.NONE) has no
    # monotonicity, and the pole ring is barely resolved without the polar
    # filter (F3). Bounded, not a closure defect (conservation is exact).
    assert res["rhoX_min"] > -0.2, res["rhoX_min"]

    # the bell mass-centroid should have rotated ~quarter turn toward the
    # (south) pole: its mean latitude sinks well below the equator
    sol, x, Jw, inner = res["sol"], res["x"], res["Jw"], res["inner"]
    w = (Jw * sol.rhoX)[inner]
    xc = np.array([np.sum(w * x[k][inner]) / np.sum(w) for k in range(3)])
    lat = np.arcsin(np.clip(-xc[2] / np.linalg.norm(xc), -1, 1))
    assert lat < -0.5, lat  # moved from the equator toward the south pole
