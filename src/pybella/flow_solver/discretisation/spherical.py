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

    #: direction of the geometric NORTH pole, fixed Cartesian components
    pole_axis_cart = (0.0, 0.0, -1.0)
    #: planetary ROTATION vector direction. NOT the north-pole direction:
    #: this embedding is a mirror image of geographic space (x2 = -r sin
    #: phi), and angular velocity is a pseudovector — eastward motion
    #: dx/dt = Omega dx/dlambda corresponds to W = +Omega z_hat in the
    #: embedded frame. Derived, and pinned by the TC2 balance gate.
    rotation_axis_cart = (0.0, 0.0, 1.0)

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

    def _e_r(self, lam, phi):
        """Radial unit vector, Cartesian components (broadcastable)."""
        cp = np.cos(phi)
        return (cp * np.cos(lam), cp * np.sin(lam), -np.sin(phi) + 0.0 * lam)

    def up_direction(self, xi):
        lam, _, phi = xi
        return self._e_r(lam, phi)

    def traditional_coriolis(self, coriolis_param):
        """``ud.coriolis_field`` callable for the traditional approximation.

        The locally-vertical component of the planetary rotation (thin
        shell / SWE): w(x) = (W . e_r) e_r with W = coriolis_param *
        rotation_axis_cart. In this mirrored embedding W . e_r =
        -coriolis_param * sin(phi) (see ``rotation_axis_cart``), so
        w_k = +c * x2 * x_k / r^2. ``coriolis_param`` uses the same
        nondimensional convention as ``ud.coriolis_strength`` (the value
        the H^-1 kernel consumes, the FULL Coriolis parameter: pass
        2*Omega_nd for a planet of rotation rate Omega). Both the sign
        and the factor are pinned empirically by the TC2 balance gate
        (the wrong sign or a factor 2 breaks geostrophy immediately).
        """
        c = float(coriolis_param)

        def field(x0, x1, x2):
            oor_sq = 1.0 / (x0**2 + x1**2 + x2**2)
            fac = c * x2 * oor_sq
            return (fac * x0, fac * x1, fac * x2)

        return field


class SphericalTerrainMap(SphericalShellMap):
    """Terrain-following radial coordinate on the sphere.

        r(lambda, eta, phi) = a + Z(eta, h(lambda, phi)),

    with Z the existing :class:`~.terrain.VerticalTransform` (Gal-Chen by
    default) on eta in [0, depth] (depth = r_top - a; the grid's radial
    coordinate axis carries eta, so ``ud.ymin = 0``, ``ud.ymax = depth``).
    Tangents by chain rule:

        t_lam = r_lam e_r + r cos(phi) e_lam,
        t_eta = r_eta e_r,
        t_phi = r_phi e_r + r e_phi,

    so the vertical coordinate lines stay RADIAL (t_eta || e_r): gravity
    remains along ``up_direction`` = e_r while the coordinate surfaces
    tilt with the terrain. For h == 0 the metric reduces to
    ``SphericalShellMap`` (bit-near; gated).

    The orography ``h(lambda, phi)`` is evaluated at lambda wrapped into
    [-pi, pi) so ghost longitudes see their periodic image (the seam
    consistency the legacy builder enforces); analytic gradient callables
    ``orography_grad = (dh_dlam, dh_dphi)`` are REQUIRED (no FD fallback:
    the grid spacing is not visible to the map). Single-component
    transforms only (SLEVE needs the smooth/residual split — later).
    """

    def __init__(self, radius, depth, orography, orography_grad, transform=None):
        super().__init__(radius, frozen_radius=False)
        self.depth = float(depth)
        if self.depth <= 0.0:
            raise ValueError("shell depth (r_top - a) must be positive")
        self.orography = orography
        self.orography_grad = orography_grad
        self.transform = (
            transform if transform is not None else terrain.GalChenTransform()
        )
        if getattr(self.transform, "n_components", 1) != 1:
            raise NotImplementedError(
                "two-component transforms (SLEVE) on the sphere need the "
                "smooth/residual orography split — not implemented yet"
            )

    @staticmethod
    def _wrap_lam(lam):
        return np.mod(lam + np.pi, 2.0 * np.pi) - np.pi

    def _r_fields(self, lam, eta, phi):
        lamw = self._wrap_lam(lam)
        h = self.orography(lamw, phi)
        dh_dlam = self.orography_grad[0](lamw, phi)
        dh_dphi = self.orography_grad[1](lamw, phi)
        tr = self.transform
        Z = tr.z(eta, h, 0.0, self.depth)
        r_eta = tr.jacobian(eta, h, 0.0, self.depth)
        r_lam = tr.slope(eta, dh_dlam, 0.0, self.depth)
        r_phi = tr.slope(eta, dh_dphi, 0.0, self.depth)
        return self.radius + Z, r_lam, r_eta, r_phi

    def coordinates(self, xi):
        lam, eta, phi = xi
        r, _, _, _ = self._r_fields(lam, eta, phi)
        cp = np.cos(phi)
        return [r * cp * np.cos(lam), r * cp * np.sin(lam), -r * np.sin(phi)]

    def tangents(self, xi):
        lam, eta, phi = xi
        r, r_lam, r_eta, r_phi = self._r_fields(lam, eta, phi)
        cl, sl = np.cos(lam), np.sin(lam)
        cp, sp = np.cos(phi), np.sin(phi)
        er = (cp * cl, cp * sl, -sp)
        elam = (-sl, cl, 0.0)
        ephi = (-sp * cl, -sp * sl, -cp)
        zero = 0.0 * (lam + eta + phi)
        t_lam = [r_lam * er[k] + r * cp * elam[k] + zero for k in range(3)]
        t_eta = [r_eta * er[k] + zero for k in range(3)]
        t_phi = [r_phi * er[k] + r * ephi[k] + zero for k in range(3)]
        return [t_lam, t_eta, t_phi]

    def height(self, xi):
        lam, eta, phi = xi
        r, _, _, _ = self._r_fields(lam, eta, phi)
        return r - self.radius
