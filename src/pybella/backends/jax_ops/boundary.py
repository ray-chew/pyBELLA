"""JAX twins of the ghost-cell / ghost-node boundary fills.

Functional, jit-compiled equivalents of
:mod:`pybella.flow_solver.utils.boundary.cell_boundary` /
``node_boundary`` / ``common.scale_wall_node_values``.

Everything structural is precomputed host-side in
:func:`build_boundary_config`: the gravity-fill index triples
(nlast/nsource/nimage) per (orientation, side, layer) — pure integer algebra
of grid shape and ghost counts — plus the per-op static factors
(stratification S at the fixed image coordinates, the dpi coefficient with
the metric-averaged dz, the ATMOSPHERIC_EXTENSION p20 differences, the
incompressible HydroState rhoY0 lookups, and oriented terrain G1/G2 slices).
The jitted kernels are then plain gathers + closed-form arithmetic.

Orientation: outside advection sweeps arrays are canonical ("phys",
vertical = ud.gravity_direction); during a sweep of the vertical axis the
arrays (and metric) are cyclically permuted so the sweep axis is last
("sweep", vertical = ndim-1). Both orientations are precomputed; terrain
arrays for "sweep" come from transposing the canonical metric — the twin
never reads the (mutated, possibly flipped) ``mem.elem.metric`` inside ops.

The gravity fill is **sequential across the two ghost layers per side**
(the outer layer's ``nlast`` is the inner ghost filled just before) — kept
as an unrolled 4-op loop, exactly like the numpy handler.
"""

import copy

import numpy as np
import jax
import jax.numpy as jnp

from pybella.utils import axes
from pybella.utils import options as opts

_MOMENTA = ("rhou", "rhov", "rhow")
_FIELDS = ("rho", "rhou", "rhov", "rhow", "rhoY", "rhoX")


# --------------------------------------------------------------------------
# no-gravity per-axis pad fills (np.pad wrap / symmetric / negative_symmetric)
# --------------------------------------------------------------------------


def _axslice(ndim, dim, sl):
    out = [slice(None)] * ndim
    out[dim] = sl
    return tuple(out)


def _pads_for(ndim, dim, ig):
    pads = [(0, 0)] * ndim
    pads[dim] = (ig, ig)
    return pads


def _inner_along(f, dim, ig):
    return f[_axslice(f.ndim, dim, slice(ig, -ig))]


def _wrap_fill(f, dim, ig):
    # jnp.pad implements the same logical modes as np.pad, exact for all
    # axis lengths (incl. the 1-inner-cell quasi-2D axis, where slice-based
    # ghost copies would read other ghosts)
    return jnp.pad(_inner_along(f, dim, ig), _pads_for(f.ndim, dim, ig), mode="wrap")


def _symmetric_fill(f, dim, ig):
    return jnp.pad(
        _inner_along(f, dim, ig), _pads_for(f.ndim, dim, ig), mode="symmetric"
    )


def _negative_symmetric_fill(f, dim, ig):
    """Exact replica of the numpy ``_negative_symmetric`` pad callback,
    operating on the zero-padded inner array with the callback's sequential
    write order (the high-side source is read AFTER the low-side write —
    they overlap on degenerate axes)."""
    v = jnp.pad(_inner_along(f, dim, ig), _pads_for(f.ndim, dim, ig))
    n = v.shape[dim]
    lo_src = -jnp.flip(v[_axslice(v.ndim, dim, slice(ig, 2 * ig))], axis=dim)
    v = v.at[_axslice(v.ndim, dim, slice(0, ig))].set(lo_src)
    hi_src = -jnp.flip(v[_axslice(v.ndim, dim, slice(n - 2 * ig, n - ig))], axis=dim)
    v = v.at[_axslice(v.ndim, dim, slice(n - ig, None))].set(hi_src)
    return v


def _no_gravity_fill(fields, dim, ig, bdry_type_int, normal_mom):
    """One axis of the no-gravity fill on the 6-field dict (traced)."""
    if bdry_type_int == _PERIODIC:
        for name in _FIELDS:
            fields[name] = _wrap_fill(fields[name], dim, ig)
    else:  # WALL
        for name in _FIELDS:
            if name == normal_mom:
                fields[name] = _negative_symmetric_fill(fields[name], dim, ig)
            else:
                fields[name] = _symmetric_fill(fields[name], dim, ig)
    return fields


_PERIODIC = 0
_WALL = 1
_GRAVITY = 2  # vertical axis with g != 0: hydrostatic fill


# --------------------------------------------------------------------------
# gravity (hydrostatic) fill: static per-op plan + jitted unrolled kernel
# --------------------------------------------------------------------------


class _GravityOp:
    """Static data for one (side, layer) hydrostatic ghost fill."""

    __slots__ = (
        "nlast",
        "nsource",
        "nimage",
        "direction",
        "S",
        "dpi_coeff",
        "dpi_static",
        "rhoY0_im",
        "G1_nsource",
        "G2_nsource",
        "G1_nimage",
        "G2_nimage",
    )


def _gravity_ops(mem, ud, y_axs, orient_perm):
    """Build the 4 sequential fill ops for one orientation.

    Replicates apply_gravity_boundary's loop structure: sides in order
    (bottom: direction=+1, offset=0), (top: direction=-1, offset=1); layers
    cur_idx in reversed(range(igs)).
    """
    elem = mem.elem
    ndim = elem.ndim
    v_phys = axes.vertical_axis(ud)
    icv = elem.sc[v_phys]
    igv = elem.igs[v_phys]
    g = ud.gravity_strength[v_phys]
    atmosphere = hasattr(ud, "ATMOSPHERIC_EXTENSION")
    metric = elem.metric

    if metric is not None:
        if not metric.vertical_line:
            raise NotImplementedError(
                "jax boundary fill for non-vertical-line (spherical) metrics "
                "is not implemented yet; run the numpy backend"
            )

        # the metric may be mid-sweep flipped when the config is first
        # built; undo the cyclic flips back to canonical (metric.vaxis
        # tracks where the physical vertical currently sits), then
        # transpose into the requested orientation
        def canonical(arr):
            a = np.asarray(arr)
            k = (v_phys - metric.vaxis) % ndim
            for _ in range(k):
                a = np.moveaxis(a, -1, 0)
            return np.transpose(a, orient_perm)

        # vertical thickness factor z_eta = J / (N_v)_v in the numpy FP
        # order (divide in the metric's own orientation, then transpose)
        z_eta = canonical(
            np.asarray(metric.J) / np.asarray(metric.N[metric.vaxis][metric.cart_v])
        )
        z = canonical(metric.height)
        G1 = canonical(metric.G1)
        G2 = canonical(metric.G2) if metric.G2 is not None else None
    else:
        y_coords = axes.coords_along(elem, v_phys)

    hydro = mem.npf.HydroState

    ops = []
    direction = -1.0
    offset = 0
    for _side in range(2):
        direction *= -1
        for cur_i in np.arange(igv)[::-1]:
            cur_idx = int(cur_i + offset * ((icv - 1) - 2 * cur_i))
            nlast = _axslice(ndim, y_axs, int(cur_idx + direction))
            nsource = _axslice(
                ndim,
                y_axs,
                int(offset * icv + direction * (2 * igv - (1 - offset) - cur_i)),
            )
            nimage = _axslice(ndim, y_axs, int(cur_idx))

            op = _GravityOp()
            op.nlast, op.nsource, op.nimage = nlast, nsource, nimage
            op.direction = float(direction)

            # stratification at the fixed image coordinates
            if metric is not None:
                op.S = jnp.asarray(1.0 / ud.stratification(z[nimage]))
            else:
                op.S = float(1.0 / ud.stratification(y_coords[nimage[y_axs]]))

            # dpi: static for ATMOSPHERIC_EXTENSION; otherwise the static
            # coefficient of (1/Y_last + S), built in the numpy FP order
            if atmosphere:
                op.dpi_static = float(
                    (hydro.p20[nimage[y_axs]] - hydro.p20[nlast[y_axs]]) * ud.Msq
                )
                op.dpi_coeff = None
            else:
                deta = elem.dxyz[v_phys]
                if metric is not None:
                    # z_eta = J / (N_v)_v, as in the numpy ghost fill
                    dz = 0.5 * (z_eta[nimage] + z_eta[nlast]) * deta
                    op.dpi_coeff = jnp.asarray(
                        direction * (mem.th.Gamma * g) * 0.5 * dz
                    )
                else:
                    op.dpi_coeff = float(direction * (mem.th.Gamma * g) * 0.5 * deta)
                op.dpi_static = None

            # incompressible branch reads the static hydrostate profile.
            # field_mode (terrain) + incompressible/atmosphere would index a
            # full field with a scalar in the numpy path — unexercised; guard.
            if hydro.field_mode and (ud.is_compressible != 1 or atmosphere):
                raise NotImplementedError(
                    "terrain hydrostates (field_mode) with incompressible or "
                    "ATMOSPHERIC_EXTENSION ghost fills are not supported on "
                    "the JAX backend (unexercised in the numpy path too)"
                )
            op.rhoY0_im = (
                float(hydro.rhoY0[nimage[y_axs]]) if ud.is_compressible != 1 else None
            )

            if metric is not None:
                op.G1_nsource = jnp.asarray(G1[nsource])
                op.G1_nimage = jnp.asarray(G1[nimage])
                op.G2_nsource = jnp.asarray(G2[nsource]) if G2 is not None else None
                op.G2_nimage = jnp.asarray(G2[nimage]) if G2 is not None else None
            else:
                op.G1_nsource = op.G1_nimage = op.G2_nsource = op.G2_nimage = None

            ops.append(op)
        offset += 1
    return ops


def _apply_gravity_ops(fields, ops, meta):
    """Unrolled sequential hydrostatic fill (traced)."""
    vert_mom = meta["vert_mom"]
    hor_moms = meta["hor_moms"]
    slope_moms = meta["slope_moms"]
    atmosphere = meta["atmosphere"]
    terrain = meta["terrain"]
    compressible = meta["compressible"]
    gm1 = meta["gm1"]
    gm1inv = meta["gm1inv"]

    for op in ops:
        rho = fields["rho"]
        rhoY = fields["rhoY"]
        rhoX = fields["rhoX"]
        vert = fields[vert_mom]
        nlast, nsource, nimage = op.nlast, op.nsource, op.nimage

        Y_last = rhoY[nlast] / rho[nlast]
        Y_source = rhoY[nsource] / rho[nsource]

        if terrain:
            slope_src = op.G1_nsource * fields[slope_moms[0]][nsource]
            if op.G2_nsource is not None:
                slope_src = slope_src + op.G2_nsource * fields[slope_moms[1]][nsource]
            contra_source = vert[nsource] - slope_src
            rhoYv_image = -contra_source * rhoY[nsource] / rho[nsource]
        else:
            rhoYv_image = -vert[nsource] * rhoY[nsource] / rho[nsource]

        S = op.S
        if atmosphere:
            dpi = op.dpi_static
        else:
            dpi = op.dpi_coeff * (1.0 / Y_last + S)

        if compressible:
            rhoY_g = ((rhoY[nlast] ** gm1) + dpi) ** gm1inv
        else:
            rhoY_g = op.rhoY0_im
        rho_g = rhoY_g * S
        Y_image = rhoY_g / rho_g

        if atmosphere:
            if op.direction > 0:  # bottom boundary
                v = vert[nsource] * Y_source / rho[nsource] * rhoY_g / Y_image
            else:  # top boundary
                v = vert[nsource] * Y_source
            Th_slc = rhoY_g / (rhoY_g / Y_image) / (rhoY[nsource] / rho[nsource])
        else:
            v = rhoYv_image / rhoY_g
            Th_slc = 1.0

        hor_vals = {m: fields[m][nsource] / rho[nsource] for m in hor_moms}
        X = rhoX[nsource] / rho[nsource]

        fields["rho"] = rho.at[nimage].set(rho_g)
        for m in hor_moms:
            fields[m] = fields[m].at[nimage].set(rho_g * hor_vals[m] * Th_slc)
        fields["rhoY"] = rhoY.at[nimage].set(rhoY_g)
        fields["rhoX"] = rhoX.at[nimage].set(rho_g * X)

        if atmosphere:
            fields[vert_mom] = fields[vert_mom].at[nimage].set(-v / (rhoY_g / rho_g))
        elif terrain:
            # slope terms read the horizontal momenta JUST assigned at nimage
            slope_im = op.G1_nimage * fields[slope_moms[0]][nimage]
            if op.G2_nimage is not None:
                slope_im = slope_im + op.G2_nimage * fields[slope_moms[1]][nimage]
            fields[vert_mom] = fields[vert_mom].at[nimage].set(rho_g * v + slope_im)
        else:
            fields[vert_mom] = fields[vert_mom].at[nimage].set(rho_g * v)

    return fields


# --------------------------------------------------------------------------
# general (non-vertical-line, spherical) metric fills
#
# Two numpy paths gain JAX twins here — both branch on
# ``metric.vertical_line is False`` and both work off the CANONICAL metric
# (`_canonical_metric`), re-oriented per orientation with the same cyclic
# rule the sweep flips use (`_orient_leaf`):
#
#   * the general free-slip WALL mirror
#     (`cell_boundary._mirror_momenta_general`): symmetric-pad the scalars,
#     then per ghost/source pair flip the wall-normal contravariant momentum
#     with the LOCAL area normals (Cramer solve, det N = J^2 > 0). Only ever
#     invoked at canonical orientation — during a sweep it is the extremal
#     (last-axis) split, where the arrays are back at 0 net flips; the
#     degenerate (thin-shell) and periodic axes never reach it mid-sweep.
#   * the well-balanced gravity ghost fill
#     (`cell_boundary._calculate_ghost_values`, e_up branch): the radial
#     free-slip wall on the compressible shell. Sequential across the two
#     ghost layers (as the vertical-line fill), phys + sweep orientations.
# --------------------------------------------------------------------------


def _canonical_metric(metric, v_phys):
    """Return the metric in canonical orientation (vaxis == v_phys).

    The config may first be built mid-sweep (the metric flipped in place);
    recover canonical on a deep copy so all orientations derive from a
    single, flip-invariant snapshot of the (time-independent) geometry."""
    ndim = np.asarray(metric.J).ndim
    k = (v_phys - metric.vaxis) % ndim
    if k == 0:
        return metric
    m = copy.deepcopy(metric)
    for _ in range(k):
        m.flip_backward()
    return m


def _orient_leaf(arr, perm):
    """Transpose a Cartesian-component leaf array into the orientation whose
    array axes are the cyclic ``perm`` of the canonical ones (matches one
    ``MetricFields.flip_forward`` composition; identity for canonical)."""
    return np.ascontiguousarray(np.transpose(np.asarray(arr), perm))


# ------------------------------------------------ general free-slip wall


def _solve_normal_system_jax(N_at, c, ndim):
    """Cramer solve of sum_k N[a][k] m_k = c_a per point (det N = J^2 > 0),
    replicating ``cell_boundary._solve_normal_system`` FP order exactly."""
    if ndim == 2:
        det = N_at[0][0] * N_at[1][1] - N_at[0][1] * N_at[1][0]
        m0 = (c[0] * N_at[1][1] - c[1] * N_at[0][1]) / det
        m1 = (N_at[0][0] * c[1] - N_at[1][0] * c[0]) / det
        return [m0, m1]
    det = (
        N_at[0][0] * (N_at[1][1] * N_at[2][2] - N_at[1][2] * N_at[2][1])
        - N_at[0][1] * (N_at[1][0] * N_at[2][2] - N_at[1][2] * N_at[2][0])
        + N_at[0][2] * (N_at[1][0] * N_at[2][1] - N_at[1][1] * N_at[2][0])
    )
    out = []
    for k in range(3):
        M = [[c[a] if kk == k else N_at[a][kk] for kk in range(3)] for a in range(3)]
        det_k = (
            M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
        )
        out.append(det_k / det)
    return out


class _GeneralWallOp:
    """Static ghost/source slices for one (side, layer) of the mirror."""

    __slots__ = ("src", "ghost")


def _general_wall_ops(shape, wall_axis, ndim, ig):
    """Ghost/source index pairs, matching ``_mirror_momenta_general``'s
    (k_layer, low/high) traversal."""
    n_cells = int(shape[wall_axis])
    ops = []
    for k_layer in range(ig):
        for low in (True, False):
            if low:
                i_ghost = k_layer
                i_src = 2 * ig - 1 - k_layer
            else:
                i_ghost = n_cells - 1 - k_layer
                i_src = n_cells - 2 * ig + k_layer
            op = _GeneralWallOp()
            op.src = _axslice(ndim, wall_axis, i_src)
            op.ghost = _axslice(ndim, wall_axis, i_ghost)
            ops.append(op)
    return ops


def _apply_general_wall(fields, N, ops, wall_axis, ig, ndim):
    """Symmetric-pad all fields (the ``_set_boundary(symmetric)`` call), then
    contravariant-mirror the momenta.

    The scalar (symmetric) pad is what the wall reflects; the momentum pad
    (tangential symmetric, rhov negsym) is fully overwritten at every ghost
    layer for a non-degenerate axis — but on the thin-shell degenerate axis
    (one interior cell) the OUTER ghost's mirror source is itself a ghost, so
    it must carry the padded value: the pad is kept, exactly as numpy does."""
    fields = _no_gravity_fill(fields, wall_axis, ig, _WALL, "rhov")
    moms = [fields[_MOMENTA[k]] for k in range(ndim)]
    for op in ops:
        sl_s, sl_g = op.src, op.ghost
        c = [
            sum(N[a][kk][sl_s] * moms[kk][sl_s] for kk in range(ndim))
            for a in range(ndim)
        ]
        c[wall_axis] = -c[wall_axis]
        N_ghost = [[N[a][kk][sl_g] for kk in range(ndim)] for a in range(ndim)]
        m_new = _solve_normal_system_jax(N_ghost, c, ndim)
        for kk in range(ndim):
            moms[kk] = moms[kk].at[sl_g].set(m_new[kk])
    for k in range(ndim):
        fields[_MOMENTA[k]] = moms[k]
    return fields


# ---------------------------------------- well-balanced gravity ghost fill


class _GeneralGravityOp:
    """Static data for one (side, layer) general (e_up) hydrostatic fill."""

    __slots__ = (
        "nlast",
        "nsource",
        "nimage",
        "direction",
        "S",
        "dpi_coeff",
        "rhoY0_im",
        "Nv_src",
        "Nv_im",
        "e_src",
        "e_im",
    )


def _general_gravity_ops(mem, ud, y_axs, perm):
    """Build the 4 sequential general (e_up) gravity ops for one orientation.

    Mirrors ``_gravity_ops``' index math; the per-op metric slices (the
    vertical area normal N_v, the up direction e_up, both at source/image)
    are static and gathered here in the requested orientation."""
    elem = mem.elem
    ndim = elem.ndim
    v_phys = axes.vertical_axis(ud)
    icv = elem.sc[v_phys]
    igv = elem.igs[v_phys]
    g = ud.gravity_strength[v_phys]
    if hasattr(ud, "ATMOSPHERIC_EXTENSION"):
        raise NotImplementedError(
            "ATMOSPHERIC_EXTENSION with a non-vertical-line (spherical) metric "
            "is not supported on the JAX backend (unexercised in numpy too)"
        )

    mc = _canonical_metric(elem.metric, v_phys)
    Nv = [_orient_leaf(mc.N[v_phys][k], perm) for k in range(ndim)]
    e_up = [_orient_leaf(mc.e_up[k], perm) for k in range(ndim)]
    height = _orient_leaf(mc.height, perm)
    h_v = _orient_leaf(mc.h_v, perm)
    hydro = mem.npf.HydroState
    if hydro.field_mode and ud.is_compressible != 1:
        raise NotImplementedError(
            "terrain hydrostates (field_mode) with an incompressible general "
            "ghost fill are not supported on the JAX backend"
        )

    ops = []
    direction = -1.0
    offset = 0
    for _side in range(2):
        direction *= -1
        for cur_i in np.arange(igv)[::-1]:
            cur_idx = int(cur_i + offset * ((icv - 1) - 2 * cur_i))
            nlast = _axslice(ndim, y_axs, int(cur_idx + direction))
            nsource = _axslice(
                ndim,
                y_axs,
                int(offset * icv + direction * (2 * igv - (1 - offset) - cur_i)),
            )
            nimage = _axslice(ndim, y_axs, int(cur_idx))

            op = _GeneralGravityOp()
            op.nlast, op.nsource, op.nimage = nlast, nsource, nimage
            op.direction = float(direction)
            op.S = jnp.asarray(1.0 / ud.stratification(height[nimage]))
            deta = elem.dxyz[v_phys]
            dz = 0.5 * (h_v[nimage] + h_v[nlast]) * deta
            op.dpi_coeff = jnp.asarray(direction * (mem.th.Gamma * g) * 0.5 * dz)
            op.rhoY0_im = (
                float(hydro.rhoY0[nimage[y_axs]]) if ud.is_compressible != 1 else None
            )
            op.Nv_src = [jnp.asarray(Nv[k][nsource]) for k in range(ndim)]
            op.Nv_im = [jnp.asarray(Nv[k][nimage]) for k in range(ndim)]
            op.e_src = [jnp.asarray(e_up[k][nsource]) for k in range(ndim)]
            op.e_im = [jnp.asarray(e_up[k][nimage]) for k in range(ndim)]
            ops.append(op)
        offset += 1
    return ops


def _apply_general_gravity_ops(fields, ops, meta):
    """Unrolled sequential general (e_up) hydrostatic fill (traced).

    Twin of ``cell_boundary._calculate_ghost_values`` / ``_assign_ghost_values``
    e_up branch: hydrostatic rho/rhoY from the along-arc pressure step, the
    tangential velocity copied per Cartesian component, and the ghost momentum
    rebuilt as (tangential part + beta e_up) with beta enforcing the
    well-balanced coordinate velocity ``v`` (rhoY flux odd across the wall)."""
    ndim = meta["ndim"]
    compressible = meta["compressible"]
    gm1 = meta["gm1"]
    gm1inv = meta["gm1inv"]

    for op in ops:
        rho = fields["rho"]
        rhoY = fields["rhoY"]
        rhoX = fields["rhoX"]
        nlast, nsource, nimage = op.nlast, op.nsource, op.nimage
        moms = [fields[_MOMENTA[k]] for k in range(ndim)]

        Y_last = rhoY[nlast] / rho[nlast]
        rho_s = rho[nsource]
        Nv_dot_m = sum(op.Nv_src[k] * moms[k][nsource] for k in range(ndim))
        u_dot_e = sum(moms[k][nsource] * op.e_src[k] for k in range(ndim)) / rho_s
        tang = [moms[k][nsource] / rho_s - u_dot_e * op.e_src[k] for k in range(ndim)]
        E_image = sum(op.Nv_im[k] * op.e_im[k] for k in range(ndim))
        Y_src = rhoY[nsource] / rho_s
        v_coord = -Y_src * Nv_dot_m / E_image

        S = op.S
        dpi = op.dpi_coeff * (1.0 / Y_last + S)
        if compressible:
            rhoY_g = (rhoY[nlast] ** gm1 + dpi) ** gm1inv
        else:
            rhoY_g = op.rhoY0_im
        rho_g = rhoY_g * S
        v = v_coord / rhoY_g
        X = rhoX[nsource] / rho_s

        fields["rho"] = rho.at[nimage].set(rho_g)
        fields["rhoY"] = rhoY.at[nimage].set(rhoY_g)
        fields["rhoX"] = rhoX.at[nimage].set(rho_g * X)
        for k in range(ndim):
            moms[k] = moms[k].at[nimage].set(rho_g * tang[k])
        Nv_dot_mt = sum(op.Nv_im[k] * moms[k][nimage] for k in range(ndim))
        beta = rho_g * v - Nv_dot_mt / E_image
        for k in range(ndim):
            moms[k] = moms[k].at[nimage].add(beta * op.e_im[k])
        for k in range(ndim):
            fields[_MOMENTA[k]] = moms[k]
    return fields


# --------------------------------------------------------------------------
# config + jitted entry kernels per (orientation, mode)
# --------------------------------------------------------------------------


class BoundaryConfig:
    """Per-(elem, ud) static boundary data + jitted fill functions."""

    def __init__(self, mem, ud):
        elem = mem.elem
        self.elem_ref = elem
        self.ud_ref = ud
        ndim = elem.ndim
        v_phys = axes.vertical_axis(ud)
        self.ndim = ndim
        self.v_phys = v_phys
        self.igs = tuple(int(i) for i in elem.igs)
        self.gravity_on = [float(ud.gravity_strength[d]) != 0.0 for d in range(ndim)]
        self.bdry_int = [
            _PERIODIC if ud.bdry_type[d] == opts.BdryType.PERIODIC else _WALL
            for d in range(ndim)
        ]
        # RAYLEIGH off the gravity axis is asserted out, matching numpy
        for d in range(ndim):
            if ud.bdry_type[d] == opts.BdryType.RAYLEIGH and not self.gravity_on[d]:
                raise AssertionError(
                    "Rayleigh boundary only defined on the gravity axis."
                )

        meta = {
            "vert_mom": _MOMENTA[v_phys],
            "hor_moms": tuple(m for i, m in enumerate(_MOMENTA) if i != v_phys),
            "atmosphere": hasattr(ud, "ATMOSPHERIC_EXTENSION"),
            "terrain": elem.metric is not None,
            "compressible": ud.is_compressible == 1,
            "gm1": mem.th.gm1,
            "gm1inv": mem.th.gm1inv,
            "slope_moms": None,
        }
        if elem.metric is not None:
            a_h1, a_h2 = axes.horizontal_axes(v_phys)
            meta["slope_moms"] = (_MOMENTA[a_h1], _MOMENTA[a_h2])
        self.meta = meta

        # general (spherical) metric: the non-vertical-line free-slip wall
        # mirror and the e_up gravity fill replace the vertical-line kernels
        self.general = elem.metric is not None and not elem.metric.vertical_line

        if any(self.gravity_on):
            # canonical ("phys") orientation: identity permutation
            ident = tuple(range(ndim))
            # sweep orientation: cyclic permutation putting v_phys last,
            # i.e. oriented axes order (v+1, ..., v) mod ndim
            sweep_perm = tuple((v_phys + 1 + i) % ndim for i in range(ndim))
            if self.general:
                phys_ops = _general_gravity_ops(mem, ud, v_phys, ident)
                sweep_ops = _general_gravity_ops(mem, ud, ndim - 1, sweep_perm)
                self.gravity_fill = {
                    "phys": self._make_general_gravity_fill(phys_ops),
                    "sweep": self._make_general_gravity_fill(sweep_ops),
                }
            else:
                phys_ops = _gravity_ops(mem, ud, v_phys, ident)
                sweep_ops = _gravity_ops(mem, ud, ndim - 1, sweep_perm)
                self.gravity_fill = {
                    "phys": self._make_gravity_fill(phys_ops),
                    "sweep": self._make_gravity_fill(sweep_ops),
                }

        self.no_gravity_fill = self._make_no_gravity_fills()

        # general free-slip WALL mirror (canonical orientation) for each
        # non-gravity WALL axis — replaces the Cartesian normal-flip fill
        self.general_wall_fill = {}
        if self.general:
            mc = _canonical_metric(elem.metric, v_phys)
            N_canon = [
                [jnp.asarray(np.asarray(mc.N[a][k])) for k in range(ndim)]
                for a in range(ndim)
            ]
            for d in range(ndim):
                if self.bdry_int[d] == _WALL and not self.gravity_on[d]:
                    self.general_wall_fill[d] = self._make_general_wall_fill(
                        N_canon, d, self.igs[d]
                    )

    def _make_general_gravity_fill(self, ops):
        meta = {
            "ndim": self.ndim,
            "compressible": self.meta["compressible"],
            "gm1": self.meta["gm1"],
            "gm1inv": self.meta["gm1inv"],
        }

        @jax.jit
        def fill(rho, rhou, rhov, rhow, rhoY, rhoX):
            fields = dict(
                rho=rho, rhou=rhou, rhov=rhov, rhow=rhow, rhoY=rhoY, rhoX=rhoX
            )
            fields = _apply_general_gravity_ops(fields, ops, meta)
            return tuple(fields[n] for n in _FIELDS)

        return fill

    def _make_general_wall_fill(self, N, wall_axis, ig):
        ndim = self.ndim

        @jax.jit
        def fill(rho, rhou, rhov, rhow, rhoY, rhoX):
            fields = dict(
                rho=rho, rhou=rhou, rhov=rhov, rhow=rhow, rhoY=rhoY, rhoX=rhoX
            )
            ops = _general_wall_ops(rho.shape, wall_axis, ndim, ig)
            fields = _apply_general_wall(fields, N, ops, wall_axis, ig, ndim)
            return tuple(fields[n] for n in _FIELDS)

        return fill

    def _make_gravity_fill(self, ops):
        meta = self.meta

        @jax.jit
        def fill(rho, rhou, rhov, rhow, rhoY, rhoX):
            fields = dict(
                rho=rho, rhou=rhou, rhov=rhov, rhow=rhow, rhoY=rhoY, rhoX=rhoX
            )
            fields = _apply_gravity_ops(fields, ops, meta)
            return tuple(fields[n] for n in _FIELDS)

        return fill

    def _make_no_gravity_fills(self):
        """One jitted single-axis fill per (dim, current_step) combination."""
        fills = {}
        for dim in range(self.ndim):
            for current_step in range(self.ndim):
                bint = self.bdry_int[current_step]
                normal_mom = _MOMENTA[current_step]
                ig = self.igs[dim]

                def make(dim=dim, bint=bint, normal_mom=normal_mom, ig=ig):
                    @jax.jit
                    def fill(rho, rhou, rhov, rhow, rhoY, rhoX):
                        fields = dict(
                            rho=rho,
                            rhou=rhou,
                            rhov=rhov,
                            rhow=rhow,
                            rhoY=rhoY,
                            rhoX=rhoX,
                        )
                        fields = _no_gravity_fill(fields, dim, ig, bint, normal_mom)
                        return tuple(fields[n] for n in _FIELDS)

                    return fill

                fills[(dim, current_step)] = make()
        return fills


_CONFIG_CACHE = {}


def get_boundary_config(mem, ud):
    key = (id(mem.elem), id(ud))
    cfg = _CONFIG_CACHE.get(key)
    if cfg is None or cfg.elem_ref is not mem.elem or cfg.ud_ref is not ud:
        cfg = BoundaryConfig(mem, ud)
        _CONFIG_CACHE[key] = cfg
    return cfg


def set_ghost_cells(mem, ud, step=None, sol=None):
    """Drop-in twin of cell_boundary.set_ghost_cells (mutates sol fields)."""
    if sol is None:
        sol = mem.sol
    cfg = get_boundary_config(mem, ud)
    ndim = cfg.ndim

    dims = range(ndim) if step is None else [ndim - 1]
    for dim in dims:
        current_step = step if step is not None else dim
        arrays = tuple(jnp.asarray(getattr(sol, n)) for n in _FIELDS)
        if cfg.gravity_on[current_step]:
            orientation = "sweep" if step is not None else "phys"
            out = cfg.gravity_fill[orientation](*arrays)
        elif dim == current_step and current_step in cfg.general_wall_fill:
            # general (spherical) free-slip wall. Only reached at canonical
            # orientation (dim == current_step): step=None fills every axis
            # canonically, and mid-sweep the only general WALL axis is the
            # extremal (last) split, where the arrays are at 0 net flips
            # (degenerate/periodic axes never reach the general mirror). The
            # mirror axis is then current_step (cell_boundary
            # ._mirror_momenta_general is called with dim=current_step).
            out = cfg.general_wall_fill[current_step](*arrays)
        else:
            out = cfg.no_gravity_fill[(dim, current_step)](*arrays)
        for name, val in zip(_FIELDS, out):
            getattr(sol, name)[...] = np.asarray(val)


# --------------------------------------------------------------------------
# node fills (set_ghost_nodes twin)
# --------------------------------------------------------------------------


def _periodic_plus_one_fill(p, dim, ig):
    """Exact replica of the numpy ``periodic_plus_one`` pad callback,
    operating on the zero-padded inner array: an overlap exchange including
    one interior node per side, both sources read before any write (np.pad
    builds the RHS tuple first). On the degenerate quasi-2D axis this
    deliberately reads transient zeros and swaps the duplicate interior
    nodes — masked by the broadcast afterwards, exactly as in numpy."""
    v = jnp.pad(_inner_along(p, dim, ig), _pads_for(p.ndim, dim, ig))
    n = v.shape[dim]
    lo_src = v[_axslice(v.ndim, dim, slice(n - 2 * ig - 1, n - ig))]
    hi_src = v[_axslice(v.ndim, dim, slice(ig, 2 * ig + 1))]
    v = v.at[_axslice(v.ndim, dim, slice(0, ig + 1))].set(lo_src)
    v = v.at[_axslice(v.ndim, dim, slice(n - ig - 1, None))].set(hi_src)
    return v


def _reflect_fill(p, dim, ig):
    return jnp.pad(_inner_along(p, dim, ig), _pads_for(p.ndim, dim, ig), mode="reflect")


_NODE_FILL_CACHE = {}


def _node_fill_fn(ndim, igs, bdry_ints, degen):
    key = (ndim, igs, bdry_ints, degen)
    fn = _NODE_FILL_CACHE.get(key)
    if fn is not None:
        return fn

    @jax.jit
    def fill(p):
        for dim in range(ndim):
            if bdry_ints[dim] == _PERIODIC:
                p = _periodic_plus_one_fill(p, dim, igs[dim])
            else:
                p = _reflect_fill(p, dim, igs[dim])
        for dim, sc in degen:
            layer = p[_axslice(p.ndim, dim, igs[dim])]
            p = jnp.broadcast_to(jnp.expand_dims(layer, axis=dim), p.shape)
        return p

    _NODE_FILL_CACHE[key] = fill
    return fill


def set_ghost_nodes(p, node, ud, igs=None):
    """Drop-in twin of node_boundary.set_ghost_nodes (mutates p)."""
    if igs is None:
        igs = node.igs
    bdry_ints = tuple(
        _PERIODIC if ud.bdry_type[d] == opts.BdryType.PERIODIC else _WALL
        for d in range(node.ndim)
    )
    degen = tuple((int(d), int(node.sc[d])) for d in axes.degenerate_axes(node))
    fn = _node_fill_fn(node.ndim, tuple(int(i) for i in igs), bdry_ints, degen)
    p[...] = np.asarray(fn(jnp.asarray(p)))


# --------------------------------------------------------------------------
# scale_wall_node_values twin
# --------------------------------------------------------------------------


def scale_wall_node_values(rhs, node, ud, factor=0.5):
    """Drop-in twin of common.scale_wall_node_values (mutates and returns rhs)."""
    ndim = node.ndim
    igs = node.igs
    out = jnp.asarray(rhs)
    for dim in range(ndim):
        is_wall = ud.bdry_type[dim] in (opts.BdryType.WALL, opts.BdryType.RAYLEIGH)
        if is_wall:
            idx = [slice(igs[d], -igs[d]) for d in range(ndim)]
            for boundary_idx in [igs[dim], -igs[dim] - 1]:
                idx[dim] = boundary_idx
                out = out.at[tuple(idx)].multiply(factor)
    rhs[...] = np.asarray(out)
    return rhs
