"""Spherical lat-lon geometry as a curvilinear map.

The sphere enters the solver the same way terrain does: as metric data
(J, area normals N_a) built from a :class:`~.terrain.CurvilinearMap`.
Momenta stay fixed GLOBAL Cartesian components; no equation, stencil or
kernel changes — only the map is new.

Computational axes (canonical orientation, matching the slab convention
h1 = x, v = y, h2 = z with ``gravity_direction = 1``):

    axis 0 (x): lambda, longitude in radians  (PERIODIC)
    axis 1 (y): r, radius (nondimensional)    (the gravity/radial axis)
    axis 2 (z): phi, latitude in radians      (WALL at +-phi_max on a
                latitude channel; POLE at +-pi/2 with ``pole=True``)

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

``frozen_radius=True`` is the thin-shell approximation: the horizontal
tangents use the fixed radius a instead of each cell's own r, so J and
every normal stop depending on r. Shallow water wants this — the layer
is thin enough that r = a everywhere in it, and the quasi-2D machinery
broadcasts one layer anyway.

It does break the identity sum_a D_a N_a = 0, which says a cell's face
vectors must cancel; without it, a uniform flow would create or destroy
mass. The leftover is exactly the radial growth term that went missing,
-2 a cos(phi) e_r, so it points straight up and the horizontal faces
still cancel perfectly. Harmless with no vertical dynamics, not harmless
with them — so the 3D compressible shell uses ``frozen_radius=False``
and keeps the identity exact. Both cases are checked in
``test_scripts/test_sphere_metric.py``.
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
    #: planetary ROTATION vector direction. Points the OPPOSITE way to
    #: the north pole above, and that is correct, not a typo: x2 = -r sin
    #: phi makes this embedding a mirror image of geographic space, and a
    #: mirror reverses the right-hand rule, so the spin arrow flips while
    #: the pole does not.
    rotation_axis_cart = (0.0, 0.0, 1.0)

    def __init__(self, radius, frozen_radius=False, pole=False):
        self.radius = float(radius)
        self.frozen_radius = bool(frozen_radius)
        self.pole = bool(pole)
        if self.radius <= 0.0:
            raise ValueError("planet radius must be positive")

    def _r_eff(self, r):
        return self.radius if self.frozen_radius else r

    def _fold_poles(self, xi):
        """Fold ghost coordinates beyond |phi| = pi/2 through the pole.

        A ghost cell past the pole covers a real point on the FAR side, at
        longitude lambda + pi. Evaluating the map at that folded point puts
        the image cell's own values in the ghost, verbatim — no rotation or
        sign flip, since momenta are global Cartesian and carry no local
        north/east basis to correct. (Negating ``t_phi`` would flip the
        sign of J, and with it every flux in that row.)

        e_r and N_phi then continue smoothly through the pole — leaning
        one way before it, upright at it, the other way after — instead of
        reflecting back on themselves as they would under a wall mirror.
        Gated by ``test_scripts/test_sphere_pole_metric.py``.

        Returns ``xi`` unchanged when ``pole`` is off (bit-identity for the
        channel cases).
        """
        if not self.pole:
            return xi
        lam, r, phi = xi[0], xi[1], xi[2]
        over = np.abs(phi) > 0.5 * np.pi
        lam_f = np.where(over, lam + np.pi, lam)
        phi_f = np.where(over, np.sign(phi) * np.pi - phi, phi)
        return [lam_f, r, phi_f]

    def pole_mask(self, xi):
        """Boolean grid mask of coordinate-singular (pole) points.

        The pole nodes |phi| = pi/2 where J = r^2 cos(phi) degenerates to
        zero; ``None`` when the map is not pole-enabled. Cell centres never
        fall exactly on the pole, so this is all-False on the cell grid.
        """
        if not self.pole:
            return None
        phi = xi[2]
        return np.abs(np.abs(phi) - 0.5 * np.pi) < 1e-9

    def coordinates(self, xi):
        lam, r, phi = self._fold_poles(xi)
        cphi = np.cos(phi)
        return [
            r * cphi * np.cos(lam),
            r * cphi * np.sin(lam),
            -r * np.sin(phi),
        ]

    def tangents(self, xi):
        lam, r, phi = self._fold_poles(xi)
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
        lam, _, phi = self._fold_poles(xi)
        return self._e_r(lam, phi)

    def traditional_coriolis(self, coriolis_param):
        """``ud.coriolis_field`` callable for the traditional approximation.

        On a thin shell only the locally-vertical part of the planetary
        rotation matters, so the rotation vector W is replaced by its
        radial projection — the usual f = 2 Omega sin(latitude):

            w(x) = (W . e_r) e_r,   W = coriolis_param * rotation_axis_cart

        ``field`` below evaluates that in Cartesian position, with no phi
        left in it, using sin(phi) = -x2/r and e_r = x/r:

            W . e_r = -c sin(phi) = c x2 / r     =>   w_k = c x2 x_k / r^2

        Pass the FULL Coriolis parameter — 2*Omega_nd for a planet
        spinning at Omega, the same convention as ``ud.coriolis_strength``
        and what the H^-1 kernel consumes. Sign and factor are both pinned
        by the Williamson-TC2 balance regression
        (``tests/test_sphere_swe_tc2.py``): either one wrong and
        geostrophy fails immediately.
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
    ``SphericalShellMap`` (bit-near; asserted in
    ``test_scripts/test_sphere_metric.py``).

    The orography ``h(lambda, phi)`` is evaluated at lambda wrapped into
    [-pi, pi) so ghost longitudes see their periodic image (the same
    periodic-seam consistency ``build_metric_fields`` enforces); analytic
    gradient callables ``orography_grad = (dh_dlam, dh_dphi)`` are
    REQUIRED (no FD fallback: the grid spacing is not visible to the
    map). Single-component transforms only (SLEVE needs the
    smooth/residual split — not implemented).
    """

    def __init__(
        self, radius, depth, orography, orography_grad, transform=None, pole=False
    ):
        super().__init__(radius, frozen_radius=False, pole=pole)
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
                "SphericalTerrainMap does not support two-component "
                "transforms yet: it passes a single orography to the "
                "transform, but SLEVE needs a (smooth, residual) pair. "
                "Use GalChenTransform."
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
        lam, eta, phi = self._fold_poles(xi)
        r, _, _, _ = self._r_fields(lam, eta, phi)
        cp = np.cos(phi)
        return [r * cp * np.cos(lam), r * cp * np.sin(lam), -r * np.sin(phi)]

    def tangents(self, xi):
        lam, eta, phi = self._fold_poles(xi)
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
        lam, eta, phi = self._fold_poles(xi)
        r, _, _, _ = self._r_fields(lam, eta, phi)
        return r - self.radius
