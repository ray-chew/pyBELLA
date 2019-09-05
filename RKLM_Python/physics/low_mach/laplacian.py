import numpy as np
from inputs.enum_bdry import BdryType
from scipy import sparse, signal
from itertools import product
from numba import jit

def stencil_9pt(elem,node,mpv,ud,diag_inv):
    igx = elem.igx
    igy = elem.igy

    icxn = node.icx
    icyn = node.icy

    iicxn = icxn - (2 * igx)
    iicyn = icyn - (2 * igy)

    iicxn, iicyn = iicyn, iicxn

    dx = node.dy
    dy = node.dx

    inner_domain = (slice(igx,-igx),slice(igy,-igy))

    hplusx = mpv.wplus[1][inner_domain].reshape(-1,)
    hplusy = mpv.wplus[0][inner_domain].reshape(-1,)

    hcenter = mpv.wcenter[inner_domain].reshape(-1,)

    diag_inv = diag_inv[inner_domain].reshape(-1,)

    oodx2 = 0.5 / (dx**2)
    oody2 = 0.5 / (dy**2)

    x_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[0] == BdryType.PERIODIC

    x_wall = ud.bdry_type[1] == BdryType.WALL
    y_wall = ud.bdry_type[0] == BdryType.WALL

    return lambda p : lap2D(p, igx,igy, iicxn, iicyn, hplusx, hplusy, hcenter, oodx2, oody2, x_periodic, y_periodic, x_wall, y_wall, diag_inv)

@jit(nopython=True, nogil=True, cache=True)
def lap2D(p, igx,igy, iicxn, iicyn, hplusx, hplusy, hcenter, oodx2, oody2, x_periodic, y_periodic, x_wall, y_wall, diag_inv):
    ngnc = (iicxn) * (iicyn)
    lap = np.zeros((ngnc))
    cnt_x = 0
    cnt_y = 0

    nine_pt = 0.25 * (2.0) * 1.0

    for idx in range(iicxn * iicyn):
        ne_topleft = idx - iicxn - 1
        ne_topright = idx - iicxn 
        ne_botleft = idx - 1
        ne_botright = idx

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

            if x_periodic:
                topmid_idx += iicxn - 1
                midmid_idx += iicxn - 1
                botmid_idx += iicxn - 1

            ne_topleft += iicxn - 1
            ne_botleft += iicxn - 1

        if cnt_x == (iicxn - 1):
            topright_idx -= iicxn - 1
            midright_idx -= iicxn - 1
            botright_idx -= iicxn - 1

            if x_periodic:
                topmid_idx -= iicxn - 1
                midmid_idx -= iicxn - 1
                botmid_idx -= iicxn - 1

            ne_topright -= iicxn - 1
            ne_botright -= iicxn - 1

        if cnt_y == 0:
            topleft_idx += ((iicxn) * (iicyn - 1)) 
            topmid_idx += ((iicxn) * (iicyn - 1))
            topright_idx += ((iicxn) * (iicyn - 1))

            if y_periodic:
                midleft_idx += ((iicxn) * (iicyn - 1))
                midmid_idx += ((iicxn) * (iicyn - 1))
                midright_idx += ((iicxn) * (iicyn - 1))

            ne_topleft += ((iicxn) * (iicyn - 1))
            ne_topright += ((iicxn) * (iicyn - 1))

        if cnt_y == (iicyn - 1):
            botleft_idx -= ((iicxn) * (iicyn - 1))
            botmid_idx -= ((iicxn) * (iicyn - 1))
            botright_idx -= ((iicxn) * (iicyn - 1))

            if y_periodic:
                midleft_idx -= ((iicxn) * (iicyn - 1))
                midmid_idx -= ((iicxn) * (iicyn - 1))
                midright_idx -= ((iicxn) * (iicyn - 1))

            ne_botleft -= ((iicxn) * (iicyn - 1))
            ne_botright -= ((iicxn) * (iicyn - 1))

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

        if x_wall and (cnt_x == 0):
            hplusx_topleft = 0.
            hplusy_topleft = 0.
            hplusx_botleft = 0.
            hplusy_botleft = 0.

        if x_wall and (cnt_x == (iicxn - 1)):
            hplusx_topright = 0.
            hplusy_topright = 0.
            hplusx_botright = 0.
            hplusy_botright = 0.

        if y_wall and (cnt_y == 0):
            hplusx_topleft = 0.
            hplusy_topleft = 0. 
            hplusx_topright = 0.
            hplusy_topright = 0.
            
        if y_wall and (cnt_y == (iicyn - 1)):
            hplusx_botleft = 0.
            hplusy_botleft = 0.  
            hplusx_botright = 0.
            hplusy_botright = 0.
                                
        dp2dxdy1 = ((midmid - midleft) - (topmid - topleft)) * nine_pt
        dp2dxdy2 = ((midright - midmid) - (topright - topmid)) * nine_pt
        dp2dxdy3 = ((botmid - botleft) - (midmid - midleft)) * nine_pt
        dp2dxdy4 = ((botright - botmid) - (midright - midmid)) * nine_pt

        lap[idx] = - hplusx_topleft * oodx2 * ((midmid - midleft) - dp2dxdy1) \
                -  hplusy_topleft * oody2 * ((midmid - topmid) - dp2dxdy1) \
                +  hplusx_topright * oodx2 * ((midright - midmid) - dp2dxdy2) \
                -  hplusy_topright * oody2 * ((midmid - topmid) + dp2dxdy2) \
                -  hplusx_botleft * oodx2 * ((midmid - midleft) + dp2dxdy3) \
                +  hplusy_botleft * oody2 * ((botmid - midmid) - dp2dxdy3) \
                +  hplusx_botright * oodx2 * ((midright - midmid) + dp2dxdy4) \
                +  hplusy_botright * oody2 * ((botmid - midmid) + dp2dxdy4) \
                +  hcenter[idx] * p[idx]

        # if cnt_x == 0 and x_wall:
        #     lap[idx] *= 2.
        # if (cnt_x == iicxn - 1) and x_wall:
        #     lap[idx] *= 2.

        # if cnt_y == 0 and y_wall:
        #     lap[idx] *= 2.
        # if (cnt_y == iicyn - 1) and y_wall:
        #     lap[idx] *= 2.

        lap[idx] *= diag_inv[idx]

        cnt_x += 1
        if cnt_x % iicxn == 0:
            cnt_y += 1
            cnt_x = 0
        
    return lap

def precon_diag_prepare(mpv, elem, node, ud):
    dx = node.dx
    dy = node.dy

    x_periodic = ud.bdry_type[0] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[1] == BdryType.PERIODIC

    igx = node.igx
    igy = node.igy

    idx_e = (slice(igx - x_periodic,-igx + x_periodic - 1),slice(igy - y_periodic, -igy + y_periodic - 1))
    idx_periodic = [slice(None)] * elem.ndim

    for dim in range(elem.ndim):
        if ud.bdry_type[dim] == BdryType.PERIODIC:
            idx_periodic[dim] = slice(1,-1)
    idx_periodic = tuple(idx_periodic)

    idx_n = (slice(igx,-igx), slice(igy,-igy))

    hplusx = mpv.wplus[0]
    hplusy = mpv.wplus[1]

    nine_pt = 0.25 * (2.0) * 1.0

    wx = (1.0 - nine_pt) / (dx**2)
    wy = (1.0 - nine_pt) / (dy**2)

    diag_kernel = np.array([[1.,1.],[1.,1.]])

    diag = np.zeros((node.sc)).squeeze()
    diag[idx_n] = -wx * signal.convolve2d(hplusx[idx_e],diag_kernel,mode='full')[idx_periodic] 
    diag[idx_n] -= wy * signal.convolve2d(hplusy[idx_e],diag_kernel,mode='full')[idx_periodic]

    diag[idx_n] += mpv.wcenter[idx_n]
    diag[idx_n] = 1.0 / diag[idx_n]

    return diag


    