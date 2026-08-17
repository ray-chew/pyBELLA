"""Hughes & Jablonowski (2023) mountain baroclinic wave — pt 2: the two
midlatitude ridges (the wave TRIGGER) via a terrain-following spherical map
+ the well-balanced "adjusted" background (Eq. 2).

pt 1 (:mod:`pybella.tests.test_hj_baroclinic`) mapped the Ullrich balanced
base state onto the SMOOTH deep shell and showed it stays (nearly) steady.
pt 2 embeds the two analytic ridges (Eq. 1, ``ullrich_baroclinic.orography``:
h0 = 2000 m at 72E / 140E, 45N) as GEOMETRY -- a ``SphericalTerrainMap``
(Gal-Chen radial coordinate) -- and re-derives the balanced state on the
tilted coordinate surfaces. No SWE bottom-topography source terms: the ridges
are geometry, and the balance is the surface-pressure adjustment (Eq. 2).

The crux is that H&J's "adjusted" balanced surface pressure (Eq. 2) is
exactly the base-state pressure profile (Eq. B4) evaluated at the surface
height z_s(lambda, phi),

    p_s(lambda, phi) = pressure(phi, z_s(lambda, phi)),

and more generally the whole well-balanced adjusted state is the analytic
Ullrich state SAMPLED at the terrain-following physical height
z = r - a = Z(eta, h(lambda, phi)) of each cell. Gal-Chen puts the bottom
coordinate surface at Z(0, h) = h = z_s, so at eta = 0 the sampled pressure
IS Eq. 2, and the vertical column above it is the analytic base state on the
lifted (tilted) surfaces. pt 1's ``sol_init`` is written ENTIRELY in
``metric.height``, so it reproduces this adjusted state unchanged once the
map carries the ridges -- the only differences from pt 1 are the map itself
and the radial axis now carrying eta in [0, depth] (so ymin = 0, ymax =
depth) instead of r.

The adjusted state is well-balanced but -- as H&J note -- NOT PERFECTLY so:
the discrete gradient-wind + hydrostatic balance on the tilted coordinate
surfaces near the ridges carries a small residual, and THAT residual is the
intended trigger. It launches the baroclinic Rossby wave train (matured over
~6 d in the paper). The coarse numpy smoke gate
(``test_scripts/test_hj_baroclinic_ridges.py``) checks only that the wave
INITIATES near the ridges (a longitude-localised meridional-wind response
well above the flat-background pt-1 adjustment) and does not blow up; the
production multi-day device run is pt 3.

As in pt 1, no initial projection: the small residual is accepted.
The now-fixed field-mode gravity ghost fill
(``cell_boundary._calculate_ghost_values``, regression
``test_scripts/test_field_mode_gravity_ghost.py``) makes the TC2-style
initial projection available should the residual prove too large -- set
``initial_projection = True``.
"""

from ..flow_solver.discretisation import spherical

from . import test_hj_baroclinic as hj
from . import ullrich_baroclinic as ub
from .case_setup import make_diag_state
from .test_hj_baroclinic import _DEPTH_M

# The IC is written entirely in ``metric.height``; with the ridge map that
# height follows the terrain, so the pt-1 initial condition IS the pt-2
# adjusted-balance state. Reuse it verbatim.
sol_init = hj.sol_init


def _nondim_orography(h_ref):
    """H&J ridge orography (Eq. 1) + analytic gradients in the map's
    NONDIMENSIONAL height units (z / h_ref).

    ``ullrich_baroclinic`` returns DIMENSIONAL metres [m] (h0 = 2000 m) and
    [m rad^-1]; a ``SphericalTerrainMap`` works in the same nondimensional
    length as its radius / depth / eta, so the surface height must be scaled
    by 1 / h_ref (the coordinates stay the angular lambda, phi in radians, so
    the gradients are only rescaled, not reparameterised). Passing the raw
    dimensional metres makes Z = eta + h b(eta) with h ~ 2000 >> depth ~ 4,
    which flips the Jacobian at the ridge latitudes.
    """
    oro = lambda lam, phi: ub.orography(lam, phi) / h_ref
    oro_grad = (
        lambda lam, phi: ub.orography_grad_lam(lam, phi) / h_ref,
        lambda lam, phi: ub.orography_grad_phi(lam, phi) / h_ref,
    )
    return oro, oro_grad


class UserData(hj.UserData):
    def __init__(self):
        super().__init__()

        a_nd = self.planet_radius / self.h_ref
        depth = _DEPTH_M / self.h_ref

        # The ridges are GEOMETRY: compose the H&J two-ridge orography (Eq. 1)
        # and its analytic gradients (required -- no FD fallback on the map)
        # into a terrain-following radial coordinate. The radial axis now
        # carries eta in [0, depth] (Gal-Chen: Z(0, h) = h, Z(depth, h) =
        # depth), not r, so ymin/ymax change with it.
        self.ymin, self.ymax = 0.0, depth
        oro, oro_grad = _nondim_orography(self.h_ref)
        self.curvilinear_map = spherical.SphericalTerrainMap(a_nd, depth, oro, oro_grad)

        self.output_base_name = "_hj_baroclinic_ridges"
        self.output_suffix = "_%i_%i" % (self.inx - 1, self.inz - 1)

        self.diag_state = make_diag_state(
            "test_hj_baroclinic_ridges",
            "target_hj_baroclinic_ridges",
            self.inx,
            self.iny,
            self.stepmax,
        )
