"""Spherical lat-lon geometry as a curvilinear map (Tier 3, sphere slice).

The sphere enters the solver the same way terrain does: as metric data
(J, area normals N_a) built from a :class:`~.terrain.CurvilinearMap`.
Momenta stay fixed GLOBAL Cartesian components; no equation, stencil or
kernel changes — only the map is new.

Computational axes (canonical orientation, matching the slab convention
h1 = x, v = y, h2 = z with ``gravity_direction = 1``):

    axis 0 (x): lambda, longitude in radians  (PERIODIC on the channel)
    axis 1 (y): r, radius (nondimensional)    (the gravity/radial axis)
    axis 2 (z): phi, latitude in radians      (WALL at +-phi_max)

Embedding into the fixed Cartesian frame (components = momenta indices):

    x_0 = r cos(phi) cos(lambda)
    x_1 = r cos(phi) sin(lambda)
    x_2 = -r sin(phi)

The north pole lies along Cartesian -x_2: with the builder's cyclic
normal convention N_a = t_{a+1} x t_{a+2} this makes the coordinate
triple (lambda, r, phi) orientation-POSITIVE, J = r~^2 cos(phi) > 0, and
N_r = J e_r (so ``e_up`` is the outward radial unit vector). The local
radial unit vector is e_r = (cos(phi)cos(lambda), cos(phi)sin(lambda),
-sin(phi)); the planetary rotation axis (toward the north pole) is the
constant Cartesian vector (0, 0, -1).

``frozen_radius=True`` evaluates the horizontal tangents at r~ = a
instead of r: J and all normals become exactly r-independent — the
thin-shell (SWE) degeneracy, consistent with the quasi-2D broadcast
machinery. The price is that the discrete metric identity sum_a D_a N_a
= 0 acquires a defect, but the defect is PURELY RADIAL (-2 a cos(phi)
e_r, from the missing d(r^2)/dr), so it never contaminates tangential
fluxes; asserted in ``test_scripts/test_sphere_metric.py``. The true
r-dependent map (``frozen_radius=False``) is what the 3D compressible
shell uses.
"""

import numpy as np

from . import terrain


class SphericalShellMap(terrain.CurvilinearMap):
    """Lat-lon-radius map of a spherical shell of reference radius ``a``.

    Parameters
    ----------
    radius : float
        Nondimensional planet radius a (ud.planet_radius / ud.h_ref).
    frozen_radius : bool
        Evaluate horizontal tangents at r~ = a (thin-shell/SWE mode; see
        module docstring). Default False (true 3D shell metric).
    """

    vertical_line = False

    #: rotation axis toward the north pole, fixed Cartesian components
    pole_axis_cart = (0.0, 0.0, -1.0)

    def __init__(self, radius, frozen_radius=False):
        self.radius = float(radius)
        self.frozen_radius = bool(frozen_radius)
        if self.radius <= 0.0:
            raise ValueError("planet radius must be positive")

    def _r_eff(self, r):
        return self.radius if self.frozen_radius else r

    def coordinates(self, xi):
        lam, r, phi = xi[0], xi[1], xi[2]
        cphi = np.cos(phi)
        return [
            r * cphi * np.cos(lam),
            r * cphi * np.sin(lam),
            -r * np.sin(phi),
        ]

    def tangents(self, xi):
        lam, r, phi = xi[0], xi[1], xi[2]
        re = self._r_eff(r)
        cl, sl = np.cos(lam), np.sin(lam)
        cp, sp = np.cos(phi), np.sin(phi)
        zero = 0.0 * (lam + r + phi)
        t_lam = [-re * cp * sl, re * cp * cl, zero]
        t_r = [cp * cl + zero, cp * sl + zero, -sp + zero]
        t_phi = [-re * sp * cl, -re * sp * sl, -re * cp + zero]
        return [t_lam, t_r, t_phi]

    def height(self, xi):
        return xi[1] - self.radius

    def traditional_coriolis(self, coriolis_param):
        """``ud.coriolis_field`` callable for the traditional approximation.

        Rotation-vector field w(x) = coriolis_param * sin(phi) * e_r —
        only the locally-vertical rotation component acts (thin shell /
        SWE). ``coriolis_param`` uses the same nondimensional convention
        as ``ud.coriolis_strength`` (the value the H^-1 kernel consumes);
        its absolute normalization (Omega vs 2*Omega) is pinned by the
        f-plane equivalence oracle, not assumed here. sin(phi) = -x2/r in
        this embedding (pole along -x2), so w_k = -c * x2 * x_k / r^2.
        """
        c = float(coriolis_param)

        def field(x0, x1, x2):
            oor_sq = 1.0 / (x0**2 + x1**2 + x2**2)
            fac = -c * x2 * oor_sq
            return (fac * x0, fac * x1, fac * x2)

        return field
