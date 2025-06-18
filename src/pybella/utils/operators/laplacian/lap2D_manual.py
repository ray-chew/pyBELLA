import numpy as np
import numba as nb
from ... import options as opts

def get_linop(mpv, node, coriolis, diag_inv, ud):
    dx = node.dx
    dy = node.dy

    hplusx = mpv.wplus[0]
    hplusy = mpv.wplus[1]
    hcenter = mpv.wcenter

    coeffs = (hplusx.T, hplusy.T, hcenter.T)

    ### Need to clean this up, but the Numba stencil is used in the Helmholtz solve for radiative BC!
    if hasattr(ud, "ATMOSPHERIC_EXTENSION") and ud.ATMOSPHERIC_EXTENSION:
        y_atmosphere = True
    else:
        y_atmosphere = False

    ###################
    x_wall = ud.bdry_type[0] == opts.BdryType.WALL
    y_wall = ud.bdry_type[1] == opts.BdryType.WALL

    cor_slc = (slice(1, -1), slice(1, -1))
    coeff_slc = (slice(1, -1), slice(1, -1))

    coeffs = (
        hplusx[coeff_slc].T.reshape(
            -1,
        ),
        hplusy[coeff_slc].T.reshape(
            -1,
        ),
        hcenter[node.i1].T.reshape(
            -1,
        ),
    )

    coriolis = (
        coriolis[0][cor_slc].reshape(
            -1,
        ),
        coriolis[1][cor_slc].reshape(
            -1,
        ),
        coriolis[2][cor_slc].reshape(
            -1,
        ),
        coriolis[3][cor_slc].reshape(
            -1,
        ),
    )

    return lambda p: lap2D_gather(
        p,
        node.iicx,
        node.iicy,
        coeffs,
        dx,
        dy,
        x_wall,
        y_wall,
        y_atmosphere,
        diag_inv[node.i1].T.reshape(
            -1,
        ),
        coriolis,
    )


@nb.njit(cache=True)
def lap2D_gather(
    p, iicxn, iicyn, coeffs, dx, dy, x_wall, y_wall, y_atmosphere, diag_inv, coriolis
):
    ngnc = (iicxn) * (iicyn)
    lap = np.zeros((ngnc))
    cnt_x = 0
    cnt_y = 0

    oodx = 1.0 / dx
    oody = 1.0 / dy
    cxx, cyy, cxy, cyx = coriolis

    hplusx, hplusy, hcenter = coeffs

    for idx in range(iicxn * iicyn):
        nr_row = idx // iicxn
        col_idx = idx - (nr_row * iicxn)

        ne_row_idx = nr_row * (iicxn + 1)
        ne_col_idx = col_idx
        ne_idx = ne_row_idx + ne_col_idx

        ne_topleft = ne_idx
        ne_topright = ne_idx + 1
        ne_botleft = ne_idx + (iicxn + 1)
        ne_botright = ne_idx + (iicxn + 1) + 1

        # get indices of the 9pt stencil
        topleft_idx = idx - iicxn - 1
        midleft_idx = idx - 1
        botleft_idx = idx + iicxn - 1

        topmid_idx = idx - iicxn
        midmid_idx = idx
        botmid_idx = idx + iicxn

        topright_idx = idx - iicxn + 1
        midright_idx = idx + 1
        botright_idx = idx + iicxn + 1

        if cnt_x == 0:
            topleft_idx += iicxn - 1
            midleft_idx += iicxn - 1
            botleft_idx += iicxn - 1

        if cnt_x == (iicxn - 1):
            topright_idx -= iicxn - 1
            midright_idx -= iicxn - 1
            botright_idx -= iicxn - 1

        val = 0

        if cnt_y == 0:
            if y_atmosphere:
                topleft_idx += 2 * (iicxn - val)
                topmid_idx += 2 * (iicxn - val)
                topright_idx += 2 * (iicxn - val)
            else:
                topleft_idx += (iicxn) * (iicyn - 1)
                topmid_idx += (iicxn) * (iicyn - 1)
                topright_idx += (iicxn) * (iicyn - 1)

        if cnt_y == (iicyn - 1):
            if y_atmosphere:
                botleft_idx -= 2 * (iicxn - val)
                botmid_idx -= 2 * (iicxn - val)
                botright_idx -= 2 * (iicxn - val)
            else:
                botleft_idx -= (iicxn) * (iicyn - 1)
                botmid_idx -= (iicxn) * (iicyn - 1)
                botright_idx -= (iicxn) * (iicyn - 1)

        topleft = p[topleft_idx]
        midleft = p[midleft_idx]
        botleft = p[botleft_idx]

        topmid = p[topmid_idx]
        midmid = p[midmid_idx]
        botmid = p[botmid_idx]

        topright = p[topright_idx]
        midright = p[midright_idx]
        botright = p[botright_idx]

        hplusx_topleft = hplusx[ne_topleft]
        hplusx_botleft = hplusx[ne_botleft]
        hplusy_topleft = hplusy[ne_topleft]
        hplusy_botleft = hplusy[ne_botleft]

        hplusx_topright = hplusx[ne_topright]
        hplusx_botright = hplusx[ne_botright]
        hplusy_topright = hplusy[ne_topright]
        hplusy_botright = hplusy[ne_botright]

        cxx_tl = cxx[ne_topleft]
        cxx_tr = cxx[ne_topright]
        cxx_bl = cxx[ne_botleft]
        cxx_br = cxx[ne_botright]

        cxy_tl = cxy[ne_topleft]
        cxy_tr = cxy[ne_topright]
        cxy_bl = cxy[ne_botleft]
        cxy_br = cxy[ne_botright]

        cyx_tl = cyx[ne_topleft]
        cyx_tr = cyx[ne_topright]
        cyx_bl = cyx[ne_botleft]
        cyx_br = cyx[ne_botright]

        cyy_tl = cyy[ne_topleft]
        cyy_tr = cyy[ne_topright]
        cyy_bl = cyy[ne_botleft]
        cyy_br = cyy[ne_botright]

        if x_wall and (cnt_x == 0):
            hplusx_topleft = 0.0
            hplusy_topleft = 0.0
            hplusx_botleft = 0.0
            hplusy_botleft = 0.0

        if x_wall and (cnt_x == (iicxn - 1)):
            hplusx_topright = 0.0
            hplusy_topright = 0.0
            hplusx_botright = 0.0
            hplusy_botright = 0.0

        if y_wall and (cnt_y == 0):
            if y_atmosphere:
                pass
            else:
                hplusx_topleft = 0.0
                hplusy_topleft = 0.0
                hplusx_topright = 0.0
                hplusy_topright = 0.0

        if y_wall and (cnt_y == (iicyn - 1)):
            if y_atmosphere:
                pass
            else:
                hplusx_botleft = 0.0
                hplusy_botleft = 0.0
                hplusx_botright = 0.0
                hplusy_botright = 0.0

        Dx_tl = 0.5 * (topmid - topleft + midmid - midleft) * hplusx_topleft
        Dx_tr = 0.5 * (topright - topmid + midright - midmid) * hplusx_topright
        Dx_bl = 0.5 * (botmid - botleft + midmid - midleft) * hplusx_botleft
        Dx_br = 0.5 * (botright - botmid + midright - midmid) * hplusx_botright

        Dy_tl = 0.5 * (midmid - topmid + midleft - topleft) * hplusy_topleft
        Dy_tr = 0.5 * (midright - topright + midmid - topmid) * hplusy_topright
        Dy_bl = 0.5 * (botmid - midmid + botleft - midleft) * hplusy_botleft
        Dy_br = 0.5 * (botright - midright + botmid - midmid) * hplusy_botright

        fac = 1.0
        Dxx = (
            0.5
            * (cxx_tr * Dx_tr - cxx_tl * Dx_tl + cxx_br * Dx_br - cxx_bl * Dx_bl)
            * oodx
            * oodx
            * fac
        )
        Dyy = (
            0.5
            * (cyy_br * Dy_br - cyy_tr * Dy_tr + cyy_bl * Dy_bl - cyy_tl * Dy_tl)
            * oody
            * oody
            * fac
        )
        Dyx = (
            0.5
            * (cxy_br * Dy_br - cxy_bl * Dy_bl + cxy_tr * Dy_tr - cxy_tl * Dy_tl)
            * oody
            * oodx
            * fac
        )
        Dxy = (
            0.5
            * (cyx_br * Dx_br - cyx_tr * Dx_tr + cyx_bl * Dx_bl - cyx_tl * Dx_tl)
            * oodx
            * oody
            * fac
        )

        lap[idx] = Dxx + Dyy + Dyx + Dxy + hcenter[idx] * p[idx]

        lap[idx] *= diag_inv[idx]

        cnt_x += 1
        if cnt_x % iicxn == 0:
            cnt_y += 1
            cnt_x = 0

    return lap
