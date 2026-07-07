"""JAX twin of :mod:`pybella.utils.operators.laplacian.lap3D`.

Same operator, with the numpy kernel's in-place steps hoisted out of the
matvec into construction time in :func:`get_linop`:

- the wall-axis coefficient slab-zeroing is data-independent and idempotent
  (the numpy kernel re-zeroes the same slabs on every call), so it is
  pre-applied to the copied coefficient boxes once;
- the periodic ghost reconstruction (tmp-swap of the duplicate rows 1 and
  n-2 plus the ghost fills) is a fixed permutation of the input along each
  periodic axis, applied in the kernel as a gather:
  ``[0 -> n-3, 1 -> n-2, n-2 -> 1, n-1 -> 2]``. (Until 2026-06-10 the
  numba ``tmp`` was a view, making the closing ``p[-2] = tmp`` a no-op;
  both backends were fixed in lockstep — proven bit-identical end-to-end,
  since the solver only visits periodically consistent vectors, on which
  mirroring and exchanging duplicate rows coincide.)

Like the numpy twin, the matvec returns the boxed (3D, padded) array — the
caller (scipy/jax bicgstab) flattens it.
"""

import numpy as np
import jax
import jax.numpy as jnp

from pybella.utils import options as opts

_COEFF = 1.0 / 16

# four node-adjacent cubes per axial direction, as in the numpy kernel
_TOPLEFTS = [
    (slice(0, None), slice(0, -1), slice(0, -1)),
    (slice(0, -1), slice(0, None), slice(0, -1)),
    (slice(0, -1), slice(0, -1), slice(0, None)),
]
_TOPRIGHTS = [
    (slice(0, None), slice(0, -1), slice(1, None)),
    (slice(1, None), slice(0, None), slice(0, -1)),
    (slice(0, -1), slice(1, None), slice(0, None)),
]
_BOTLEFTS = [
    (slice(0, None), slice(1, None), slice(0, -1)),
    (slice(0, -1), slice(0, None), slice(1, None)),
    (slice(1, None), slice(0, -1), slice(0, None)),
]
_BOTRIGHTS = [
    (slice(0, None), slice(1, None), slice(1, None)),
    (slice(1, None), slice(0, None), slice(1, None)),
    (slice(1, None), slice(1, None), slice(0, None)),
]


def _four_cube_sum(q, axis):
    return (
        q[_TOPLEFTS[axis]]
        + q[_TOPRIGHTS[axis]]
        + q[_BOTLEFTS[axis]]
        + q[_BOTRIGHTS[axis]]
    )


def _difference_block(c, flx, axis):
    """D_axis(c * flx): node-summed minus/plus parts of the differenced flux."""
    q = c * flx
    sl_m = [slice(None)] * 3
    sl_p = [slice(None)] * 3
    sl_m[axis] = slice(None, -1)
    sl_p[axis] = slice(1, None)
    qm = _four_cube_sum(q[tuple(sl_m)], axis)
    qp = _four_cube_sum(q[tuple(sl_p)], axis)
    return qm, qp


def get_linop(elem, node, npf, ud, diag_inv, dt, cij):
    """Build the full-tensor 27-point operator matvec on the node.isc box.

    Same contract as the numpy twin; see its docstring for the
    discretisation. The returned callable maps the C-order ravel of the
    node.isc box to the boxed result array.
    """
    oodxyz = 1.0 / (node.dxyz**2)
    oodx2, oody2, oodz2 = oodxyz[0], oodxyz[1], oodxyz[2]
    odx, ody, odz = 1.0 / node.dx, 1.0 / node.dy, 1.0 / node.dz

    i1 = (slice(1, -1), slice(1, -1), slice(1, -1))

    ndim = elem.ndim
    periodic = [ud.bdry_type[dim] == opts.BdryType.PERIODIC for dim in range(ndim)]

    C = [[np.ascontiguousarray(cij[i][j][i1]) for j in range(3)] for i in range(3)]

    # cross blocks only enter when H^-1 has off-diagonal content; decided on
    # the un-zeroed boxes, exactly as the numpy twin does at construction
    use_cross = bool(
        max(np.max(np.abs(C[i][j])) for i in range(3) for j in range(3) if i != j) > 0.0
    )

    # pre-apply the wall-axis coefficient slab-zeroing the numpy kernel
    # performs (idempotently) on every call
    for dim in range(ndim):
        if not periodic[dim]:
            lo = [slice(None)] * 3
            hi = [slice(None)] * 3
            lo[dim] = 0
            hi[dim] = -1
            for i in range(3):
                for j in range(3):
                    C[i][j][tuple(lo)] = 0.0
                    C[i][j][tuple(hi)] = 0.0

    hcenter = np.ascontiguousarray(npf.wcenter[i1])
    diag_inv = np.ascontiguousarray(diag_inv)

    shx, shy, shz = hcenter.shape
    padded = (shx + 2, shy + 2, shz + 2)

    # periodic ghost reconstruction as a fixed permutation per axis
    perms = []
    for dim in range(3):
        n = padded[dim]
        perm = np.arange(n)
        if dim < ndim and periodic[dim]:
            perm[[0, 1, -2, -1]] = [n - 3, n - 2, 1, 2]
        perms.append(perm)

    # tree_util.Partial, not a closure (see lap2D.get_linop): `padded` is
    # re-derived from hcenter's static shape inside the trace, and use_cross
    # selects between two module-level variants, so both stay static
    apply_fn = _lap3D_apply_cross if use_cross else _lap3D_apply_plain
    return jax.tree_util.Partial(
        apply_fn,
        C=tuple(tuple(jnp.asarray(C[i][j]) for j in range(3)) for i in range(3)),
        hcenter=jnp.asarray(hcenter),
        scales=(oodx2, oody2, oodz2, odx, ody, odz),
        perms=tuple(jnp.asarray(p) for p in perms),
        diag_inv=jnp.asarray(diag_inv),
    )


# --------------------------------------------------------------------------
# pole-ring collapse (Stage F, F7): Galerkin master-node embedding, twin of
# ``numerics.pole_collapse.PoleCollapse``. The wrapped operator
# ``gather o A o scatter`` solves each pole ring as a single master unknown;
# the non-master ring entries stay exactly zero (scatter ignores them, gather
# re-zeroes them), so the bicgstab plumbing is untouched.
# --------------------------------------------------------------------------


def _collapsed_apply(v, raw, scatter_src, ring_complement, uniq_mem, uniq_master):
    """gather(raw(scatter(v))) on the flat solve vector (traced).

    The ring is zeroed by a MASK MULTIPLY (``* ring_complement``), not a
    ``.at[ring_all].set(0.0)`` integer scatter — jax's bicgstab wraps the
    operator in ``custom_linear_solve``, which double-transposes it, and an
    integer scatter-SET is not double-transposable (a non-unique scatter-ADD
    is). Same numerics as ``numerics.pole_collapse.PoleCollapse.gather``."""
    v = jnp.asarray(v, dtype=jnp.float64)
    y = jnp.reshape(raw(v[scatter_src]), (-1,))
    contrib = y[uniq_mem]  # read BEFORE masking (uniq_mem lands in the ring)
    y = y * ring_complement
    y = y.at[uniq_master].add(contrib)
    return y


def wrap_pole_collapse(raw, coll):
    """Wrap a lap3D matvec (a ``tree_util.Partial``) with the pole-ring
    Galerkin scatter/gather from a numpy ``PoleCollapse`` (index arrays only).

    Returns a ``tree_util.Partial`` so repeated solves keep a stable function
    identity and hit ``elliptic_solve._solve``'s jit cache (see its docstring
    on the OOM risk of per-call closures)."""
    ring_complement = np.ones(coll.n)
    ring_complement[coll.ring_all] = 0.0
    return jax.tree_util.Partial(
        _collapsed_apply,
        raw=raw,
        scatter_src=jnp.asarray(coll.scatter_src),
        ring_complement=jnp.asarray(ring_complement),
        uniq_mem=jnp.asarray(coll.uniq_mem),
        uniq_master=jnp.asarray(coll.uniq_master),
    )


def _lap3D_apply(p, C, hcenter, scales, perms, diag_inv, use_cross):
    # jnp.asarray (not np.asarray): must accept jax tracers as well as
    # eager numpy probes from the equivalence tests
    p = jnp.asarray(p, dtype=jnp.float64)
    padded = tuple(s + 2 for s in hcenter.shape)
    return _lap3D(p, C, hcenter, scales, padded, perms, diag_inv, use_cross)


@jax.jit
def _lap3D_apply_plain(p, C, hcenter, scales, perms, diag_inv):
    return _lap3D_apply(p, C, hcenter, scales, perms, diag_inv, False)


@jax.jit
def _lap3D_apply_cross(p, C, hcenter, scales, perms, diag_inv):
    return _lap3D_apply(p, C, hcenter, scales, perms, diag_inv, True)


def _lap3D(p0, C, hcenter, scales, padded, perms, diag_inv, use_cross):
    oodx2, oody2, oodz2, odx, ody, odz = scales

    p = p0.reshape(padded)
    p = jnp.take(p, perms[0], axis=0)
    p = jnp.take(p, perms[1], axis=1)
    p = jnp.take(p, perms[2], axis=2)

    # cell-averaged directional differences F_j(p) on the cell box
    x_fluxes = p[1:, :, :] - p[:-1, :, :]
    y_fluxes = p[:, 1:, :] - p[:, :-1, :]
    z_fluxes = p[:, :, 1:] - p[:, :, :-1]

    x_flx = _four_cube_sum(x_fluxes, 0)
    y_flx = _four_cube_sum(y_fluxes, 1)
    z_flx = _four_cube_sum(z_fluxes, 2)

    # diagonal blocks: D_i(C_ii F_i)
    x_flxm, x_flxp = _difference_block(C[0][0], x_flx, 0)
    y_flxm, y_flxp = _difference_block(C[1][1], y_flx, 1)
    z_flxm, z_flxp = _difference_block(C[2][2], z_flx, 2)

    interior = (
        oodx2 * _COEFF * (-x_flxm + x_flxp)
        + oody2 * _COEFF * (-y_flxm + y_flxp)
        + oodz2 * _COEFF * (-z_flxm + z_flxp)
        + hcenter * p[1:-1, 1:-1, 1:-1]
    )

    if use_cross:
        ods = (odx, ody, odz)
        flxs = (x_flx, y_flx, z_flx)
        cross = jnp.zeros_like(x_flxm)
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                qm, qp = _difference_block(C[i][j], flxs[j], i)
                cross = cross + ods[i] * ods[j] * _COEFF * (qp - qm)
        interior = interior + cross

    lap = jnp.zeros(padded, dtype=p.dtype).at[1:-1, 1:-1, 1:-1].set(interior)

    return lap * diag_inv
