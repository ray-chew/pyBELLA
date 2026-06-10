"""Oracle B: the full-H^-1 3D elliptic path vs the 2D path, WITH Coriolis.

A y-uniform 3D x-z problem with rotation Omega_3 = (Ox, Oy, Oz) is
physically identical to the 2D x-y problem obtained by relabeling z -> y.
That relabeling is an odd axis permutation, so the rotation PSEUDOvector
maps with a sign flip:

    Omega_2D = -(Ox, Oz, Oy),    (u, v, w)_2D = (u, w, v)_3D.

Since the 2D elliptic path embeds the full H^-1 coefficients in its stencil
(trusted, regression-locked), this validates the 3D full-tensor operator's
cross terms — signs, axis pairing, and layout — against it. Requires g = 0
and nonhydrostasy = 1 (the hydrostatic mask and buoyancy term sit on
different slots under the odd swap).

Uses the travelling-vortex setup of test_3d_elliptic_oracle with
Omega_3 = (100, 0, 100) (the 3d-coriolis regression case's rotation,
omega * dt = 1: strong coupling).
"""

import numpy as np

from pybella.flow_solver.numerics import implicit_euler

from test_3d_elliptic_oracle import build_ud, build_state

TOL = 1e-6
OMEGA = 100.0


def test_3d_full_coriolis_matches_2d():
    ud3 = build_ud(iny=2, inz=65)
    ud3.coriolis_strength = np.array([OMEGA, 0.0, OMEGA])
    mem3 = build_state(ud3, "xz")
    implicit_euler.do_implicit_part(mem3, ud3, dt=0.01)

    ud2 = build_ud(iny=65, inz=1)
    # odd permutation (z -> y): pseudovector sign flip
    ud2.coriolis_strength = np.array([-OMEGA, -OMEGA, 0.0])
    mem2 = build_state(ud2, "xy")
    implicit_euler.do_implicit_part(mem2, ud2, dt=0.01)

    jc, jn = 2, 2  # interior y cell / node slice of the 3D state
    checks = {
        "p2_nodes": (mem3.npf.p2_nodes[:, jn, :], mem2.npf.p2_nodes),
        "rhou": (mem3.sol.rhou[:, jc, :], mem2.sol.rhou),
        "rhow vs rhov (vertical)": (mem3.sol.rhow[:, jc, :], mem2.sol.rhov),
        "rhov vs rhow (out-of-plane)": (mem3.sol.rhov[:, jc, :], mem2.sol.rhow),
        "p2 y-uniformity": (mem3.npf.p2_nodes[:, 2, :], mem3.npf.p2_nodes[:, 3, :]),
        "rhou y-uniformity": (mem3.sol.rhou[:, 2, :], mem3.sol.rhou[:, 3, :]),
    }

    failures = []
    for name, (lhs, rhs) in checks.items():
        diff = np.max(np.abs(lhs - rhs))
        if diff >= TOL:
            failures.append(f"{name}: max|diff| = {diff:.3e}")
    assert not failures, "; ".join(failures)
