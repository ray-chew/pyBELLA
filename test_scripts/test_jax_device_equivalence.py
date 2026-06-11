"""Device-resident step equivalence (phase B of the device plan).

Validation ladder below the full-run gate (test_jax_device_fullrun.py):

1. multi-step window equivalence — ``device_step.run_window`` vs the hybrid
   ``time_update.do`` on deep-copied states, all sol fields + p2_nodes at
   magnitude-scaled 1e-5 (the established Krylov convergence-slack floor;
   the two paths solve ulp-different elliptic systems);
2. compile audit — exactly 2 step compilations (Strang parity) over 3 steps;
3. seam test — ``time_update.do`` itself routes to the device driver under
   ``ud.backend = "jax-device"``;
4. guard test — unsupported configs raise with an actionable message.

Skips cleanly when jax is not installed.
"""

import copy

import numpy as np
import pytest

jax = pytest.importorskip("jax")

from pybella.backends.jax_ops import device_step  # noqa: E402
from pybella.flow_solver.discretisation import time_update  # noqa: E402

import jax_equiv_fixtures as fx  # noqa: E402

TOL = 1e-5  # Krylov convergence-slack floor (see component-2 measurement)


class _NullDebug:
    def write(self, *a, **k):
        pass

    def populate(self, *a, **k):
        pass


def assert_close(got, want, tol=TOL, label=""):
    want = np.asarray(want)
    scale = max(1.0, float(np.max(np.abs(want))))
    diff = float(np.max(np.abs(np.asarray(got) - want)))
    assert (
        diff <= tol * scale
    ), f"{label}: max|diff| = {diff:.3e} > {tol:.0e} * {scale:.3e}"


def _window_case(mem, ud, nsteps=3):
    ud.aux = getattr(ud, "aux", "") or ""
    ud.stepmax = nsteps
    tout = 1e9  # run exactly nsteps steps

    mem_h = copy.deepcopy(mem)
    mem_d = copy.deepcopy(mem)
    try:
        ud.backend = "jax"
        time_update.do(mem_h, ud, tout, None, None, _NullDebug())
        ud.backend = "jax-device"
        device_step.run_window(mem_d, ud, tout, writer=None)
    finally:
        ud.backend = "numpy"

    assert mem_d.time.step == mem_h.time.step == nsteps
    assert abs(mem_d.time.t - mem_h.time.t) < 1e-14
    for name in ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX"):
        assert_close(getattr(mem_d.sol, name), getattr(mem_h.sol, name), label=name)
    assert_close(mem_d.npf.p2_nodes, mem_h.npf.p2_nodes, label="p2_nodes")
    return mem_d


def test_window_igw_wall():
    mem, ud = fx.make_igw_mem()
    mem_d = _window_case(mem, ud)
    # compile audit: parity 0 and 1 only
    assert mem_d._device_compile_count == 2


def test_window_lamb_atmosphere_rayleigh():
    mem, ud = fx.make_lamb_mem()
    _window_case(mem, ud)


def test_window_unstable_lamb_forcing():
    from pybella.tests import test_unstable_lamb

    mem, ud = fx._mem_from_case(test_unstable_lamb)
    assert ud.rayleigh_forcing
    _window_case(mem, ud, nsteps=2)


def test_window_agnesi3d_terrain():
    mem, ud = fx.make_agnesi3d_mem(True)
    _window_case(mem, ud, nsteps=2)


def test_window_agnesi2d_terrain():
    mem, ud = fx.make_agnesi2d_mem()
    _window_case(mem, ud, nsteps=2)


def test_seam_routes_to_device():
    from unittest import mock

    mem, ud = fx.make_igw_mem()
    ud.stepmax = 1
    called = {"n": 0}
    orig = device_step.run_window

    def spy(*a, **k):
        called["n"] += 1
        return orig(*a, **k)

    try:
        ud.backend = "jax-device"
        with mock.patch.object(device_step, "run_window", side_effect=spy):
            time_update.do(mem, ud, 1e9, None, None, _NullDebug())
    finally:
        ud.backend = "numpy"
    assert called["n"] == 1


def test_guard_raises_on_blending():
    mem, ud = fx.make_igw_mem()
    ud.continuous_blending = True
    try:
        ud.backend = "jax-device"
        with pytest.raises(NotImplementedError, match="blending"):
            device_step.run_window(mem, ud, 1e9, writer=None)
    finally:
        ud.backend = "numpy"
        ud.continuous_blending = False
