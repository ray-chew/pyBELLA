"""Small shared helpers for regression-case ``UserData`` construction.

The physics in each ``tests/test_*.py`` ``UserData`` is case-specific and stays
in its own module. These pure helpers factor out only the two genuinely-repeated,
error-prone mechanical idioms so they are written once:

* ``build_bdry`` — the per-instance, object-dtype boundary-type triple (must be a
  fresh array per instance, never a shared class attribute).
* ``make_diag_state`` — the ``DiagnosticState`` wiring, which centralises the
  ``Nx = inx - 1`` / ``Ny = iny - 1`` / ``steps = [stepmax - 1]`` offsets that are
  easy to mistype, while forwarding any case-specific keywords (``plot_compare``,
  ``time_increment``, ``tolerances`` ...).

Neither helper mutates global state; each returns a value the case assigns.
"""

from types import SimpleNamespace

import numpy as np

from ..utils import axes
from ..utils import options as opts
from ..utils.data_structures import DiagnosticState


def build_bdry(x_bc, y_bc, z_bc):
    """Return a fresh ``(3,)`` object-dtype array of boundary types."""
    bdry = np.empty((3), dtype=object)
    bdry[0] = x_bc
    bdry[1] = y_bc
    bdry[2] = z_bc
    return bdry


def make_diag_state(test_name, file_name, inx, iny, stepmax, steps=None, **kwargs):
    """Build a ``DiagnosticState`` with the standard index/step offsets."""
    return DiagnosticState(
        test_name=test_name,
        file_name=file_name,
        Nx=inx - 1,
        Ny=iny - 1,
        steps=steps if steps is not None else [stepmax - 1],
        **kwargs,
    )


def mirror_centers(coord, c, cm):
    """Per-coordinate nearest periodic image of a vortex centre.

    For each entry of ``coord`` pick whichever of the centre ``c`` or its
    domain-wrapped image ``cm`` is closer, so the radius used by a vortex IC is
    measured against the nearest image on a periodic domain. Returns a fresh
    array shaped like ``coord``.
    """
    cc = np.zeros_like(coord)
    cc[...] = c * (np.abs(coord - c) < np.abs(coord - cm))
    cc[...] += cm * (np.abs(coord - c) > np.abs(coord - cm))
    return cc


def apply_rayleigh_bdry(ud, elem=None, node=None, *, with_tau=False):
    """Switch the vertical-axis boundary to RAYLEIGH when ``rayleigh_bdry_switch``.

    The vertical axis comes from ``axes.vertical_axis(ud)``, never a hard-coded
    index (axial-agnosticity invariant). With ``with_tau`` the sponge profiles
    ``ud.tcy, ud.tny`` are also set in the same block via
    ``rayleigh_boundary.get_tau_y`` (the terrain cases); the Lamb cases set their
    sponge separately and pass ``with_tau=False``. A no-op unless
    ``ud.rayleigh_bdry_switch`` is set.
    """
    if not getattr(ud, "rayleigh_bdry_switch", False):
        return
    ud.bdry_type[axes.vertical_axis(ud)] = opts.BdryType.RAYLEIGH
    if with_tau:
        from ..flow_solver.utils.boundary import rayleigh_boundary

        ud.tcy, ud.tny = rayleigh_boundary.get_tau_y(ud, elem, node, 0.5)


def do_initial_projection(Sol, npf, elem, node, th, ud, *, u0, v0, w0=0.0):
    """Project a balanced vortex IC onto the discrete incompressible constraint.

    Freeze the regime to incompressible, subtract the background wind
    ``(u0, v0, w0)``, run one implicit-Euler step, then restore the nodal
    pressure and re-add the wind.
    A no-op unless ``ud.initial_projection``. For 2D cases ``w0 == 0.0`` makes
    the ``rhow`` terms an exact (``±0``) no-op, so the same helper serves 2D and
    3D.
    """
    if ud.initial_projection != True:
        return

    from ..flow_solver.numerics import implicit_euler
    from ..flow_solver.utils import cache

    is_compressible = np.copy(ud.is_compressible)
    compressibility = np.copy(ud.compressibility)
    ud.is_compressible = 0
    ud.compressibility = 0.0

    p2aux = np.copy(npf.p2_nodes)

    Sol.rhou -= u0 * Sol.rho
    Sol.rhov -= v0 * Sol.rho
    Sol.rhow -= w0 * Sol.rho

    mem = SimpleNamespace()
    mem.sol = Sol
    mem.npf = npf
    mem.elem = elem
    mem.node = node
    mem.th = th
    mem.time = SimpleNamespace()
    mem.time.t = ud.dtfixed
    mem.time.step = 0
    mem.cache = cache.FlowSolverCache()

    implicit_euler.do_implicit_part(
        mem, ud, ud.dtfixed, writer=None, label="initial_projection"
    )

    npf.p2_nodes[...] = p2aux
    npf.dp2_nodes[...] = 0.0

    Sol.rhou += u0 * Sol.rho
    Sol.rhov += v0 * Sol.rho
    Sol.rhow += w0 * Sol.rho

    ud.is_compressible = is_compressible
    ud.compressibility = compressibility
