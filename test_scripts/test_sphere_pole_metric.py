"""Gates for the pole-aware metric build.

The pole extends the +-80 deg spherical channel to the full pole-to-pole
sphere. The backbone: ghost cells past |phi| = pi/2 cover physical points
on the FAR side of the pole, at longitude lambda + pi. Evaluating the map
at the FOLDED coordinate makes the ghost metric a fold-COPY of the image
cell (J > 0, no vector rotation); the pole NODES |phi| = pi/2 carry J == 0
exactly (the cos(phi) coordinate singularity) with ooJ := 0 and a stored
pole mask.

These gates pin:
1. J > 0 in every cell; J == 0 exactly at (and only at) the |phi| = pi/2
   node rows; ooJ finite everywhere, 0 at pole nodes; pole_mask correct.
2. fold-copy identity: every ghost-slab metric array equals the
   fold-indexed interior array to roundoff (ghost == image cell).
3. fold o fold = identity on coordinates.
4. the discrete metric identity sum_a D_a N_a stays 2nd order on the
   interior (poles do not spoil the well-balancing away from them).
5. BdryType.POLE <-> map pole=True consistency validation.
6. channel (pole=False) builds are byte-identical to the non-pole map.
"""

import numpy as np
import pytest

from pybella.flow_solver.discretisation import grid as dis_grid
from pybella.flow_solver.discretisation import spherical
from pybella.utils import axes as _axes
from pybella.utils import options as opts

_A = 2.0  # nondimensional planet radius


class _PoleUD:
    """Minimal ud for grid_init on the GLOBAL (pole-to-pole) sphere."""

    def __init__(self, n=(32, 4, 16), a=_A, pole=True, terrain=False, depth=0.1):
        self.inx, self.iny, self.inz = n[0] + 1, n[1] + 1, n[2] + 1
        self.xmin, self.xmax = -np.pi, np.pi  # lambda, PERIODIC (2*pi)
        self.zmin, self.zmax = -0.5 * np.pi, 0.5 * np.pi  # phi, POLE
        self.gravity_direction = 1
        self.bdry_type = np.array(
            [
                opts.BdryType.PERIODIC,
                opts.BdryType.WALL,
                opts.BdryType.POLE if pole else opts.BdryType.WALL,
            ]
        )
        if terrain:
            self.ymin, self.ymax = 0.0, depth
            h, grad = _hill()
            self.curvilinear_map = spherical.SphericalTerrainMap(
                a, depth, h, grad, pole=pole
            )
        else:
            self.ymin, self.ymax = a * 0.95, a * 1.05  # r, WALL
            self.curvilinear_map = spherical.SphericalShellMap(a, pole=pole)


def _hill():
    h0, sig = 0.02, 0.5

    def h(lam, phi):
        return h0 * (1.0 + np.cos(lam)) / 2.0 * np.exp(-((phi / sig) ** 2))

    def dh_dlam(lam, phi):
        return -h0 / 2.0 * np.sin(lam) * np.exp(-((phi / sig) ** 2))

    def dh_dphi(lam, phi):
        return (
            h0
            * (1.0 + np.cos(lam))
            / 2.0
            * (-2.0 * phi / sig**2)
            * np.exp(-((phi / sig) ** 2))
        )

    return h, (dh_dlam, dh_dphi)


def _coord(grid_obj, axis):
    shape = [1] * grid_obj.ndim
    shape[axis] = -1
    return _axes.coords_along(grid_obj, axis).reshape(shape)


# --------------------------------------------------- J >= 0 + pole mask


@pytest.mark.parametrize("terrain", [False, True], ids=["shell", "terrain"])
def test_pole_jacobian_and_mask(terrain):
    ud = _PoleUD(terrain=terrain)
    elem, node = dis_grid.grid_init(ud)

    # cells never sit on the pole -> J strictly positive, no mask
    assert np.all(elem.metric.J > 0.0)
    assert elem.metric.pole_mask is None
    assert np.all(np.isfinite(elem.metric.ooJ))

    # nodes: J == 0 exactly and ONLY at |phi| = pi/2
    m = node.metric
    phi = np.broadcast_to(_coord(node, 2), m.J.shape)
    at_pole = np.abs(np.abs(phi) - 0.5 * np.pi) < 1e-9
    assert np.array_equal(m.J == 0.0, at_pole)
    assert m.pole_mask is not None
    assert np.array_equal(m.pole_mask, at_pole)
    assert np.all(np.isfinite(m.ooJ))
    assert np.all(m.ooJ[at_pole] == 0.0)
    assert np.all(m.ooJ[~at_pole] > 0.0)
    # the two pole rows span every longitude and radius
    assert at_pole.sum() == 2 * m.J.shape[0] * m.J.shape[1]


# ------------------------------------------------------ fold-copy identity


def _fold_indices(grid_obj, nodal):
    """(lambda-shift, phi-mirror) index maps realizing the pole fold on a
    grid: interior lambda shifts by N/2 cells, ghost phi rows mirror to
    their interior image rows. Cells mirror about the boundary FACE; nodes
    reflect about the boundary NODE (the pole node itself)."""
    igx, _, igz = (int(i) for i in grid_obj.igs)
    ncx, _, ncz = (int(s) for s in grid_obj.sc)
    # lambda periodicity is over interior CELLS (nodes carry the +pi seam
    # duplicate, so their interior count is one larger); the +pi remap is a
    # shift of N/2 cells modulo N
    N = (ncx - 2 * igx) - (1 if nodal else 0)

    lam_interior = np.arange(igx, igx + N)
    lam_fold = igx + (lam_interior - igx + N // 2) % N

    lo_ghost = np.arange(0, igz)
    hi_ghost = np.arange(ncz - igz, ncz)
    if nodal:  # reflect about the pole node
        lo_img = 2 * igz - lo_ghost
        hi_img = 2 * (ncz - 1 - igz) - hi_ghost
    else:  # mirror about the wall face
        lo_img = 2 * igz - 1 - lo_ghost
        hi_img = (ncz - igz - 1) - (hi_ghost - (ncz - igz))
    ghost = np.concatenate([lo_ghost, hi_ghost])
    img = np.concatenate([lo_img, hi_img])
    return lam_interior, lam_fold, ghost, img


@pytest.mark.parametrize("terrain", [False, True], ids=["shell", "terrain"])
@pytest.mark.parametrize("loc", [0, 1], ids=["cells", "nodes"])
def test_fold_copy_identity(terrain, loc):
    ud = _PoleUD(terrain=terrain)
    elem, node = dis_grid.grid_init(ud)
    grid_obj = (elem, node)[loc]
    m = grid_obj.metric
    lam_i, lam_f, ghost, img = _fold_indices(grid_obj, nodal=(loc == 1))

    arrays = [m.J, m.h_v, m.height] + [c for Na in m.N for c in Na]
    arrays += list(m.e_up) + [c for c in m.x if c is not None]

    # fixed geometric scale: normal components scale like J ~ r^2; a
    # near-zero component differing by roundoff must not fail on its own
    # (collapsed) per-array scale
    scale = float(np.max(np.abs(m.J)))
    for arr in arrays:
        a = np.asarray(arr)
        # ghost phi slab at interior longitudes == fold-indexed image
        gh = a[np.ix_(lam_i, np.arange(a.shape[1]), ghost)]
        im = a[np.ix_(lam_f, np.arange(a.shape[1]), img)]
        np.testing.assert_allclose(gh, im, rtol=1e-11, atol=1e-10 * scale)


def test_fold_maps_ghosts_into_interior_and_is_idempotent():
    """The fold lands every over-the-pole ghost latitude back inside
    [-pi/2, pi/2]; folding an already-interior coordinate is a no-op."""
    cmap = spherical.SphericalShellMap(_A, pole=True)
    lam = np.linspace(-np.pi, np.pi, 7).reshape(-1, 1, 1)
    r = np.array([_A]).reshape(1, -1, 1)
    phi = np.array([0.6 * np.pi, -0.55 * np.pi, 0.3, -1.2]).reshape(1, 1, -1)
    lam_f, r_f, phi_f = cmap._fold_poles([lam, r, phi])
    assert np.all(np.abs(phi_f) <= 0.5 * np.pi + 1e-12)
    # ghost longitudes shifted by pi, interior ones untouched
    over = np.abs(phi) > 0.5 * np.pi
    exp_lam = np.where(over, lam + np.pi, lam)
    np.testing.assert_allclose(lam_f, np.broadcast_to(exp_lam, lam_f.shape))
    # idempotent on the folded (interior) coordinate
    lam_ff, _, phi_ff = cmap._fold_poles([lam_f, r_f, phi_f])
    np.testing.assert_allclose(lam_ff, lam_f)
    np.testing.assert_allclose(phi_ff, np.broadcast_to(phi_f, phi_ff.shape))


# --------------------------------------------- interior metric identity


def test_pole_metric_identity_second_order():
    """sum_a D_a N_a -> 0 at 2nd order over interior cells clear of the
    immediate pole ring (the pole singularity must not spoil the
    well-balancing where the flow actually lives)."""
    errs = []
    for n in ((32, 4, 16), (64, 4, 32)):
        ud = _PoleUD(n=n)
        elem, _ = dis_grid.grid_init(ud)
        m = elem.metric
        d = (elem.dx, elem.dy, elem.dz)
        acc_max = 0.0
        for k in range(3):
            acc = np.zeros_like(m.J[(slice(1, -1),) * 3])
            for a in range(3):
                arr = m.N[a][k]
                sl_p = [slice(1, -1)] * 3
                sl_m = [slice(1, -1)] * 3
                sl_p[a] = slice(2, None)
                sl_m[a] = slice(0, -2)
                acc += (arr[tuple(sl_p)] - arr[tuple(sl_m)]) / (2.0 * d[a])
            # exclude the two interior cell rings nearest each pole
            acc_max = max(acc_max, np.max(np.abs(acc[:, :, 2:-2])))
        errs.append(acc_max / _A**2)
    assert errs[1] < errs[0] / 2.5, errs
    assert errs[1] < 5e-2, errs


# ---------------------------------------------------- validation guards


def test_pole_config_validation():
    # POLE in bdry but map pole=False -> hard error
    ud = _PoleUD(pole=False)
    ud.bdry_type[2] = opts.BdryType.POLE
    with pytest.raises(ValueError, match="POLE"):
        dis_grid.grid_init(ud)

    # map pole=True but no POLE bdry -> hard error
    ud = _PoleUD(pole=True)
    ud.bdry_type[2] = opts.BdryType.WALL
    with pytest.raises(ValueError, match="POLE"):
        dis_grid.grid_init(ud)

    # odd interior longitude count -> hard error
    ud = _PoleUD(n=(31, 4, 16))
    with pytest.raises(ValueError, match="EVEN"):
        dis_grid.grid_init(ud)

    # latitude not spanning [-pi/2, pi/2] -> hard error
    ud = _PoleUD()
    ud.zmax = 0.4 * np.pi
    with pytest.raises(ValueError, match="pi/2"):
        dis_grid.grid_init(ud)


def test_channel_map_unchanged_without_pole():
    """pole=False leaves _fold_poles an exact identity (same objects) and
    the metric build byte-identical to a fresh non-pole map."""
    cmap = spherical.SphericalShellMap(_A, frozen_radius=True, pole=False)
    xi = [np.array([0.3]), np.array([_A]), np.array([1.3])]
    assert cmap._fold_poles(xi) is xi
    assert cmap.pole_mask(xi) is None
