"""JAX twin of :mod:`pybella.utils.operators.laplacian.lap2D_manual`.

The numpy kernel loops over flat node indices, computing 9-point-stencil
neighbour indices with periodic / atmosphere wrapping and zeroing wall
coefficients per node. All of that depends only on the grid shape and the
static BC flags — never on the data — so :func:`get_linop` precomputes:

- the nine gather-index arrays of the stencil (x-corrections applied before
  y-corrections, additively, exactly as in the scalar code);
- the per-corner coefficient and Coriolis vectors, gathered once and with
  the wall zeroing pre-applied as boolean masks.

The jitted matvec is then pure vectorized arithmetic over nine gathers,
term-for-term identical to the scalar loop body.
"""

import numpy as np
import jax
import jax.numpy as jnp

from pybella.utils import options as opts


def get_linop(npf, node, coriolis, diag_inv, ud):
    """2D (x-y plane) stencil operator; same contract as the numpy twin."""
    dx = node.dx
    dy = node.dy

    y_atmosphere = bool(
        hasattr(ud, "ATMOSPHERIC_EXTENSION") and ud.ATMOSPHERIC_EXTENSION
    )

    # RAYLEIGH is a sponged wall, exactly as in the numpy twin
    _wall = (opts.BdryType.WALL, opts.BdryType.RAYLEIGH)
    x_wall = ud.bdry_type[0] in _wall
    y_wall = ud.bdry_type[1] in _wall

    cor_slc = (slice(1, -1), slice(1, -1))
    coeff_slc = (slice(1, -1), slice(1, -1))

    hplusx = np.ravel(npf.wplus[0][coeff_slc], order="F")
    hplusy = np.ravel(npf.wplus[1][coeff_slc], order="F")
    hcenter = np.ravel(npf.wcenter[node.i1], order="F")

    cxx = np.ravel(coriolis[0][cor_slc], order="C")
    cyy = np.ravel(coriolis[1][cor_slc], order="C")
    cxy = np.ravel(coriolis[2][cor_slc], order="C")
    cyx = np.ravel(coriolis[3][cor_slc], order="C")

    dinv = np.ravel(diag_inv[node.i1], order="F")

    iicxn, iicyn = node.iicx, node.iicy
    N = iicxn * iicyn
    idx = np.arange(N)
    cnt_y, cnt_x = np.divmod(idx, iicxn)

    # --- 9-pt stencil gather indices, with the scalar loop's wrap
    # corrections applied in the same order (x first, then y, additive) ---
    stencil = {
        "topleft": idx - iicxn - 1,
        "midleft": idx - 1,
        "botleft": idx + iicxn - 1,
        "topmid": idx - iicxn,
        "midmid": idx.copy(),
        "botmid": idx + iicxn,
        "topright": idx - iicxn + 1,
        "midright": idx + 1,
        "botright": idx + iicxn + 1,
    }

    left = cnt_x == 0
    right = cnt_x == iicxn - 1
    top = cnt_y == 0
    bot = cnt_y == iicyn - 1

    for name in ("topleft", "midleft", "botleft"):
        stencil[name][left] += iicxn - 1
    for name in ("topright", "midright", "botright"):
        stencil[name][right] -= iicxn - 1

    # (the scalar code's `val` is always 0: top/bottom wrap is 2*iicxn rows
    # in atmosphere mode, a full grid otherwise)
    y_wrap = 2 * iicxn if y_atmosphere else iicxn * (iicyn - 1)
    for name in ("topleft", "topmid", "topright"):
        stencil[name][top] += y_wrap
    for name in ("botleft", "botmid", "botright"):
        stencil[name][bot] -= y_wrap

    # --- per-corner node->edge coefficient gathers, wall-masked ---
    ne_idx = cnt_y * (iicxn + 1) + cnt_x
    ne = {
        "tl": ne_idx,
        "tr": ne_idx + 1,
        "bl": ne_idx + (iicxn + 1),
        "br": ne_idx + (iicxn + 1) + 1,
    }

    hpx = {c: hplusx[ne[c]] for c in ne}
    hpy = {c: hplusy[ne[c]] for c in ne}

    if x_wall:
        for c in ("tl", "bl"):
            hpx[c][left] = 0.0
            hpy[c][left] = 0.0
        for c in ("tr", "br"):
            hpx[c][right] = 0.0
            hpy[c][right] = 0.0
    if y_wall and not y_atmosphere:
        for c in ("tl", "tr"):
            hpx[c][top] = 0.0
            hpy[c][top] = 0.0
        for c in ("bl", "br"):
            hpx[c][bot] = 0.0
            hpy[c][bot] = 0.0

    cor = {
        name: {c: arr[ne[c]] for c in ne}
        for name, arr in (("cxx", cxx), ("cyy", cyy), ("cxy", cxy), ("cyx", cyx))
    }

    # tree_util.Partial, not a closure: the stable function identity plus
    # coefficient-arrays-as-pytree-leaves let elliptic_solve._solve reuse
    # one compiled bicgstab across steps (see its docstring)
    return jax.tree_util.Partial(
        _lap2D_apply,
        stencil={k: jnp.asarray(v) for k, v in stencil.items()},
        hpx={k: jnp.asarray(v) for k, v in hpx.items()},
        hpy={k: jnp.asarray(v) for k, v in hpy.items()},
        cor={k: {c: jnp.asarray(v) for c, v in d.items()} for k, d in cor.items()},
        hcenter=jnp.asarray(hcenter),
        dinv=jnp.asarray(dinv),
        oodx=1.0 / dx,
        oody=1.0 / dy,
    )


@jax.jit
def _lap2D_apply(p, stencil, hpx, hpy, cor, hcenter, dinv, oodx, oody):
    # jnp.asarray (not np.asarray): must accept jax tracers as well as
    # eager numpy probes from the equivalence tests
    p = jnp.asarray(p, dtype=jnp.float64)
    return _lap2D_gather(p, stencil, hpx, hpy, cor, hcenter, dinv, oodx, oody)


def _lap2D_gather(p, stencil, hpx, hpy, cor, hcenter, dinv, oodx, oody):
    topleft = p[stencil["topleft"]]
    midleft = p[stencil["midleft"]]
    botleft = p[stencil["botleft"]]
    topmid = p[stencil["topmid"]]
    midmid = p[stencil["midmid"]]
    botmid = p[stencil["botmid"]]
    topright = p[stencil["topright"]]
    midright = p[stencil["midright"]]
    botright = p[stencil["botright"]]

    Dx_tl = 0.5 * (topmid - topleft + midmid - midleft) * hpx["tl"]
    Dx_tr = 0.5 * (topright - topmid + midright - midmid) * hpx["tr"]
    Dx_bl = 0.5 * (botmid - botleft + midmid - midleft) * hpx["bl"]
    Dx_br = 0.5 * (botright - botmid + midright - midmid) * hpx["br"]

    Dy_tl = 0.5 * (midmid - topmid + midleft - topleft) * hpy["tl"]
    Dy_tr = 0.5 * (midright - topright + midmid - topmid) * hpy["tr"]
    Dy_bl = 0.5 * (botmid - midmid + botleft - midleft) * hpy["bl"]
    Dy_br = 0.5 * (botright - midright + botmid - midmid) * hpy["br"]

    cxx, cyy, cxy, cyx = cor["cxx"], cor["cyy"], cor["cxy"], cor["cyx"]

    fac = 1.0
    Dxx = (
        0.5
        * (
            cxx["tr"] * Dx_tr
            - cxx["tl"] * Dx_tl
            + cxx["br"] * Dx_br
            - cxx["bl"] * Dx_bl
        )
        * oodx
        * oodx
        * fac
    )
    Dyy = (
        0.5
        * (
            cyy["br"] * Dy_br
            - cyy["tr"] * Dy_tr
            + cyy["bl"] * Dy_bl
            - cyy["tl"] * Dy_tl
        )
        * oody
        * oody
        * fac
    )
    Dyx = (
        0.5
        * (
            cxy["br"] * Dy_br
            - cxy["bl"] * Dy_bl
            + cxy["tr"] * Dy_tr
            - cxy["tl"] * Dy_tl
        )
        * oody
        * oodx
        * fac
    )
    Dxy = (
        0.5
        * (
            cyx["br"] * Dx_br
            - cyx["tr"] * Dx_tr
            + cyx["bl"] * Dx_bl
            - cyx["tl"] * Dx_tl
        )
        * oodx
        * oody
        * fac
    )

    return (Dxx + Dyy + Dyx + Dxy + hcenter * midmid) * dinv
