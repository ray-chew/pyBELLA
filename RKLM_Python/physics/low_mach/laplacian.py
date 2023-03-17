import numpy as np
from inputs.enum_bdry import BdryType
from scipy import sparse, signal
from numba import jit, njit, prange
import numba as nb

def stencil_9pt(elem,node,mpv,Sol,ud,diag_inv,dt,coriolis_params):
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
    i1 = node.i1

    hplusx = mpv.wplus[1][i1].reshape(-1,)
    hplusy = mpv.wplus[0][i1].reshape(-1,)
    hcenter = mpv.wcenter[i1].reshape(-1,)

    diag_inv = diag_inv[i1].reshape(-1,)

    oodx = 1.0 / (dx)
    oody = 1.0 / (dy)

    x_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[0] == BdryType.PERIODIC

    x_wall = ud.bdry_type[1] == BdryType.WALL or ud.bdry_type[1] == BdryType.RAYLEIGH
    y_wall = ud.bdry_type[0] == BdryType.WALL or ud.bdry_type[0] == BdryType.RAYLEIGH


    #### Compute Coriolis parameters:
    # nonhydro = ud.nonhydrostasy
    # g = ud.gravity_strength[1]
    # Msq = ud.Msq
    # dy = elem.dy

    # wh1, wv, wh2 = dt * ud.coriolis_strength

    # igs = elem.igs

    # ndim = node.ndim

    # nindim = np.empty((ndim),dtype='object')
    # innerdim = np.copy(nindim)
    # innerdim1 = np.copy(nindim)
    # eindim = np.empty((ndim),dtype='object')

    # for dim in range(ndim):
    #     is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
    #     nindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
    #     innerdim[dim] = slice(igs[dim],-igs[dim])
    #     eindim[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)

    #     if dim == 1:
    #         y_idx = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)
    #         right_idx = None if -igs[dim]+is_periodic == 0 else -igs[dim]+is_periodic
    #         y_idx1 = slice(igs[dim]-is_periodic+1, right_idx)

    #     innerdim1[dim] = slice(igs[dim]-1, (-igs[dim]+1))
 
    # strat = (mpv.HydroState_n.S0[y_idx1] - mpv.HydroState_n.S0[y_idx]) / dy

    # nindim = tuple(nindim)
    # eindim = tuple(eindim)
    # innerdim = tuple(innerdim)
    # innerdim1 = tuple(innerdim1)

    # for dim in range(0,elem.ndim,2):
    #     is_periodic = ud.bdry_type[dim] != BdryType.PERIODIC
    #     strat = np.expand_dims(strat, dim)
    #     strat = np.repeat(strat, elem.sc[dim]-int(2*is_periodic+igs[dim]), axis=dim)

    # Y = Sol.rhoY[nindim] / Sol.rho[nindim]

    # nu = np.zeros_like(mpv.wcenter)
    # nu[eindim] = -dt**2 * (g / Msq) * strat * Y

    # nu = nu[inner_domain].reshape(-1,)

    # denom = 1.0 / (wh1**2 + wh2**2 + (nu + nonhydro) * (wv**2 + 1))

    # coeff_uu = (wh1**2 + nu + nonhydro) * denom
    # coeff_uv = nonhydro * (wh1 * wv + wh2) * denom
    # coeff_vu = (wh1 * wv - wh2) * denom
    # coeff_vv = nonhydro * (1 + wv**2) * denom

    # coriolis_params = (coeff_vv, coeff_vu, coeff_uv, coeff_uu)
    # coriolis_params = np.array(coriolis_params, dtype=np.float64)

    return lambda p : lap2D_gather(p, igx,igy, iicxn, iicyn, hplusx, hplusy, hcenter, oodx, oody, x_periodic, y_periodic, x_wall, y_wall, diag_inv, coriolis_params)

    # ndim = elem.ndim
    # periodicity = np.zeros(ndim, dtype='int')
    # for dim in range(0,ndim,2):
    #     periodicity[dim] = ud.bdry_type[dim] == BdryType.PERIODIC

    # # proj = (slice(None,),slice(None,),igz)
    # i0 = (slice(0,-1),slice(0,-1))
    # i1 = (slice(1,-1),slice(1,-1))
    # i2 = (slice(igx,-igx),slice(igy,-igy))

    # hplusx = mpv.wplus[0][i0][i1]
    # hplusy = mpv.wplus[1][i0][i1]
    # hcenter = mpv.wcenter[i2]
    # diag_inv = diag_inv[i1]

    # dxy = (dx, dy)

    # #### Compute Coriolis parameters:
    # nonhydro = ud.nonhydrostasy
    # g = ud.gravity_strength[1]
    # Msq = ud.Msq
    # dy = elem.dy

    # wh1, wv, wh2 = dt * ud.coriolis_strength

    # first_nodes_row_right_idx = (slice(1,None))
    # first_nodes_row_left_idx = (slice(0,-1))

    # strat = (mpv.HydroState_n.S0[first_nodes_row_right_idx] - mpv.HydroState_n.S0[first_nodes_row_left_idx]) / dy

    # for dim in range(0,elem.ndim,2):
    #     strat = np.expand_dims(strat, dim)
    #     strat = np.repeat(strat, elem.sc[dim], axis=dim)

    # Y = Sol.rhoY / Sol.rho
    # nu = -dt**2 * (g / Msq) * strat * Y

    # coeff_uu = (wh1**2 + nu + nonhydro)
    # coeff_uv = nonhydro * (wh1 * wv + wh2)
    # coeff_vu = (wh1 * wv - wh2)
    # coeff_vv = nonhydro * (1 + wv**2)

    # coriolis_params = (coeff_uu[i1], coeff_uv, coeff_vu, coeff_vv)

    # return lambda p : lap2Dc(p, hplusx, hplusy, hcenter, dxy, periodicity, diag_inv, coriolis_params)

@jit(nopython=True, nogil=False, cache=True)
def lap2D_gather(p, igx,igy, iicxn, iicyn, hplusx, hplusy, hcenter, oodx, oody, x_periodic, y_periodic, x_wall, y_wall, diag_inv, coriolis):
    ngnc = (iicxn) * (iicyn)
    lap = np.zeros((ngnc))
    cnt_x = 0
    cnt_y = 0

    nine_pt = 0.25 * (2.0) * 1.0
    cyy, cxx, cyx, cxy = coriolis
    oodx2 = 0.5 * oodx**2
    oody2 = 0.5 * oody**2

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

            ne_topleft += iicxn - 1
            ne_botleft += iicxn - 1

        if cnt_x == (iicxn - 1):
            topright_idx -= iicxn - 1
            midright_idx -= iicxn - 1
            botright_idx -= iicxn - 1

            ne_topright -= iicxn - 1
            ne_botright -= iicxn - 1

        if cnt_y == 0:
            topleft_idx += ((iicxn) * (iicyn - 1)) 
            topmid_idx += ((iicxn) * (iicyn - 1))
            topright_idx += ((iicxn) * (iicyn - 1))

            ne_topleft += ((iicxn) * (iicyn - 1))
            ne_topright += ((iicxn) * (iicyn - 1))

        if cnt_y == (iicyn - 1):
            botleft_idx -= ((iicxn) * (iicyn - 1))
            botmid_idx -= ((iicxn) * (iicyn - 1))
            botright_idx -= ((iicxn) * (iicyn - 1))

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

        cxx_tl  = cxx[ne_topleft]
        cxx_tr = cxx[ne_topright]
        cxx_bl  = cxx[ne_botleft]
        cxx_br = cxx[ne_botright]

        cxy_tl  = cxy[ne_topleft]
        cxy_tr  = cxy[ne_topright]
        cxy_bl  = cxy[ne_botleft]
        cxy_br  = cxy[ne_botright]

        cyx_tl  = cyx[ne_topleft]
        cyx_tr  = cyx[ne_topright]
        cyx_bl  = cyx[ne_botleft]
        cyx_br  = cyx[ne_botright]

        cyy_tl  = cyy[ne_topleft]
        cyy_tr  = cyy[ne_topright]
        cyy_bl  = cyy[ne_botleft]
        cyy_br  = cyy[ne_botright]
        

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
                                
        # dp2dxdy1 = ((midmid - midleft) - (topmid - topleft)) * nine_pt
        # dp2dxdy2 = ((midright - midmid) - (topright - topmid)) * nine_pt
        # dp2dxdy3 = ((botmid - botleft) - (midmid - midleft)) * nine_pt
        # dp2dxdy4 = ((botright - botmid) - (midright - midmid)) * nine_pt

        # lap[idx] = - hplusx_topleft * oodx2 * ((midmid - midleft) - dp2dxdy1) \
        #         -  hplusy_topleft * oody2 * ((midmid - topmid) - dp2dxdy1) \
        #         +  hplusx_topright * oodx2 * ((midright - midmid) - dp2dxdy2) \
        #         -  hplusy_topright * oody2 * ((midmid - topmid) + dp2dxdy2) \
        #         -  hplusx_botleft * oodx2 * ((midmid - midleft) + dp2dxdy3) \
        #         +  hplusy_botleft * oody2 * ((botmid - midmid) - dp2dxdy3) \
        #         +  hplusx_botright * oodx2 * ((midright - midmid) + dp2dxdy4) \
        #         +  hplusy_botright * oody2 * ((botmid - midmid) + dp2dxdy4) \
        #         +  hcenter[idx] * p[idx]

        Dx_tl = 0.5 * (topmid   - topleft + midmid   - midleft) * hplusx_topleft
        Dx_tr = 0.5 * (topright - topmid  + midright - midmid ) * hplusx_topright
        Dx_bl = 0.5 * (botmid   - botleft + midmid   - midleft) * hplusx_botleft
        Dx_br = 0.5 * (botright - botmid  + midright - midmid ) * hplusx_botright

        Dy_tl = 0.5 * (midmid   - topmid   + midleft - topleft) * hplusy_topleft
        Dy_tr = 0.5 * (midright - topright + midmid  - topmid ) * hplusy_topright
        Dy_bl = 0.5 * (botmid   - midmid   + botleft - midleft) * hplusy_botleft
        Dy_br = 0.5 * (botright - midright + botmid  - midmid ) * hplusy_botright

        fac = 1.0
        Dxx = 0.5 * (cxx_tr * Dx_tr - cxx_tl * Dx_tl + cxx_br * Dx_br - cxx_bl * Dx_bl) * oodx * oodx * fac
        Dyy = 0.5 * (cyy_br * Dy_br - cyy_tr * Dy_tr + cyy_bl * Dy_bl - cyy_tl * Dy_tl) * oody * oody * fac
        Dyx = 0.5 * (cyx_br * Dy_br - cyx_bl * Dy_bl + cyx_tr * Dy_tr - cyx_tl * Dy_tl) * oody * oodx * fac
        Dxy = 0.5 * (cxy_br * Dx_br - cxy_tr * Dx_tr + cxy_bl * Dy_bl - cxy_tl * Dx_tl) * oodx * oody * fac
        

        lap[idx] = Dxx + Dyy + Dyx + Dxy + hcenter[idx] * p[idx]


        # fac = 1.0
        # lap[idx] = - fac * hplusx_topleft * oodx2 * ((midmid - midleft)) \
        #         -  fac * hplusy_topleft * oody2 * ((midmid - topmid)) \
        #         +  fac * hplusx_topright * oodx2 * ((midright - midmid)) \
        #         -  fac * hplusy_topright * oody2 * ((midmid - topmid)) \
        #         -  fac * hplusx_botleft * oodx2 * ((midmid - midleft)) \
        #         +  fac * hplusy_botleft * oody2 * ((botmid - midmid)) \
        #         +  fac * hplusx_botright * oodx2 * ((midright - midmid)) \
        #         +  fac * hplusy_botright * oody2 * ((botmid - midmid)) \
        #         +  hcenter[idx] * p[idx]

        lap[idx] *= diag_inv[idx]

        cnt_x += 1
        if cnt_x % iicxn == 0:
            cnt_y += 1
            cnt_x = 0
        
    return lap


@jit(nopython=True, nogil=False, cache=True)
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

            # if x_periodic:
            #     topmid_idx += iicxn - 1
            #     midmid_idx += iicxn - 1
            #     botmid_idx += iicxn - 1

            ne_topleft += iicxn - 1
            ne_botleft += iicxn - 1

        if cnt_x == (iicxn - 1):
            topright_idx -= iicxn - 1
            midright_idx -= iicxn - 1
            botright_idx -= iicxn - 1

            # if x_periodic:
            #     topmid_idx -= iicxn - 1
            #     midmid_idx -= iicxn - 1
            #     botmid_idx -= iicxn - 1

            ne_topright -= iicxn - 1
            ne_botright -= iicxn - 1

        if cnt_y == 0:
            topleft_idx += ((iicxn) * (iicyn - 1)) 
            topmid_idx += ((iicxn) * (iicyn - 1))
            topright_idx += ((iicxn) * (iicyn - 1))

            # if y_periodic:
            #     midleft_idx += ((iicxn) * (iicyn - 1))
            #     midmid_idx += ((iicxn) * (iicyn - 1))
            #     midright_idx += ((iicxn) * (iicyn - 1))

            ne_topleft += ((iicxn) * (iicyn - 1))
            ne_topright += ((iicxn) * (iicyn - 1))

        if cnt_y == (iicyn - 1):
            botleft_idx -= ((iicxn) * (iicyn - 1))
            botmid_idx -= ((iicxn) * (iicyn - 1))
            botright_idx -= ((iicxn) * (iicyn - 1))

            # if y_periodic:
            #     midleft_idx -= ((iicxn) * (iicyn - 1))
            #     midmid_idx -= ((iicxn) * (iicyn - 1))
            #     midright_idx -= ((iicxn) * (iicyn - 1))

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


def stencil_9pt_numba_test(mpv,node,coriolis,diag_inv, ud):
    dx = node.dx
    dy = node.dy

    hplusx = mpv.wplus[0]
    hplusy = mpv.wplus[1]
    hcenter = mpv.wcenter

    coeffs = (hplusx.T, hplusy.T, hcenter.T)

    shp = node.iisc

    dummy_p = np.zeros((node.isc[1],node.isc[0]))

    # coeff_slc = (slice(1,None), slice(1,None))
    # coeffs = (hplusx[coeff_slc].T, hplusy[coeff_slc].T, hcenter.T)
    # cor_slc = (slice(1,None), slice(1,None))
    # coriolis = (coriolis[0][cor_slc],coriolis[1][cor_slc],coriolis[2][cor_slc],coriolis[3][cor_slc])

    if hasattr(ud, 'LAMB_BDRY'):
        return lambda p : lap2D_numba_test(p, dummy_p, dx, dy, coeffs, diag_inv.T, coriolis, shp)

    ###################
    else:
        x_wall = ud.bdry_type[0] == BdryType.WALL or ud.bdry_type[0] == BdryType.RAYLEIGH
        y_wall = ud.bdry_type[1] == BdryType.WALL or ud.bdry_type[1] == BdryType.RAYLEIGH

        y_rayleigh = ud.bdry_type[1] == BdryType.RAYLEIGH

        cor_slc = (slice(1,-1), slice(1,-1))
        # cor_slc = (slice(0,-2), slice(0,-2))
        coeff_slc = (slice(1,-1), slice(1,-1))

        # cor_slc = (slice(None,), slice(None,))
        # coeff_slc = (slice(None,), slice(None,))

        coeffs = (hplusx[coeff_slc].T.reshape(-1,), hplusy[coeff_slc].T.reshape(-1,), hcenter[node.i1].T.reshape(-1,))
        
        coriolis = (coriolis[0][cor_slc].reshape(-1,),coriolis[1][cor_slc].reshape(-1,),coriolis[2][cor_slc].reshape(-1,),coriolis[3][cor_slc].reshape(-1,))

        return lambda p : lap2D_gather_new(p, node.iicx, node.iicy, coeffs, dx, dy, y_rayleigh, x_wall, y_wall, diag_inv[node.i1].T.reshape(-1,), coriolis)
    

@jit(nopython=True, cache=False)
def lap2D_numba_test(p, dp, dx, dy, coeffs, diag_inv, coriolis, shp):
    p = p.reshape(shp[1],shp[0])
    dp[1:-1,1:-1] = p

    # dp = p
    dp = periodic(dp)
    dp = kernel_9pt(dp, dx, dy, coeffs[0], coeffs[1], coeffs[2], diag_inv, coriolis[0], coriolis[1], coriolis[2], coriolis[3])
    # p = dp
    # dp = periodic(dp)
    p = dp[1:-1,1:-1]

    return p.ravel()


@njit(cache=True)
def lap2D_gather_new(p, iicxn, iicyn, coeffs, dx, dy, y_rayleigh, x_wall, y_wall, diag_inv, coriolis):
    ngnc = (iicxn) * (iicyn)
    lap = np.zeros((ngnc))
    cnt_x = 0
    cnt_y = 0

    oodx = 1.0 / dx
    oody = 1.0 / dy
    cxx, cyy, cxy, cyx = coriolis

    hplusx, hplusy, hcenter = coeffs

    for idx in range(iicxn * iicyn):
        # ne_topleft = idx - (iicxn - 1)
        # ne_topright = idx - (iicxn)
        # ne_botleft = idx - 1
        # ne_botright = idx

        nr_row = idx // iicxn
        col_idx = idx - (nr_row * iicxn)

        ne_row_idx = nr_row * (iicxn + 1)
        ne_col_idx = col_idx
        ne_idx = ne_row_idx + ne_col_idx

        # ne_topleft = ne_idx - (iicxn + 1) - 1
        # ne_topright = ne_idx - (iicxn + 1)
        # ne_botleft = ne_idx - 1
        # ne_botright = ne_idx

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

            # ne_topleft += iicxn - 1
            # ne_botleft += iicxn - 1

            # ne_topleft += iicxn + 1
            # ne_botleft += iicxn + 1

        if cnt_x == (iicxn - 1):
            topright_idx -= iicxn - 1
            midright_idx -= iicxn - 1
            botright_idx -= iicxn - 1

            # ne_topright -= iicxn - 1
            # ne_botright -= iicxn - 1

            # ne_topright -= iicxn + 1
            # ne_botright -= iicxn + 1

        if cnt_y == 0:
            topleft_idx += ((iicxn) * (iicyn - 1)) 
            topmid_idx += ((iicxn) * (iicyn - 1))
            topright_idx += ((iicxn) * (iicyn - 1))

            # ne_topleft += ((iicxn + 1) * (iicyn))
            # ne_topright += ((iicxn + 1) * (iicyn))

            # ne_topleft += ((iicxn + 2) * (iicyn + 1))
            # ne_topright += ((iicxn + 2) * (iicyn + 1))

        # if cnt_y == (iicyn - 1) and not y_rayleigh:
        if cnt_y == (iicyn - 1):
            botleft_idx -= ((iicxn) * (iicyn - 1))
            botmid_idx -= ((iicxn) * (iicyn - 1))
            botright_idx -= ((iicxn) * (iicyn - 1))

            # ne_botleft -= ((iicxn + 1) * (iicyn))
            # ne_botright -= ((iicxn + 1) * (iicyn))

            # ne_botleft -= ((iicxn + 2) * (iicyn + 1))
            # ne_botright -= ((iicxn + 2) * (iicyn + 1))

        ############
        # top BC handling for rayleigh
        ############
        # if cnt_y == (iicyn - 1) and y_rayleigh:
        #     botleft_idx = topleft_idx
        #     botmid_idx = topmid_idx
        #     botright_idx = topright_idx

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

        cxx_tl  = cxx[ne_topleft]
        cxx_tr = cxx[ne_topright]
        cxx_bl  = cxx[ne_botleft]
        cxx_br = cxx[ne_botright]

        cxy_tl  = cxy[ne_topleft]
        cxy_tr  = cxy[ne_topright]
        cxy_bl  = cxy[ne_botleft]
        cxy_br  = cxy[ne_botright]

        cyx_tl  = cyx[ne_topleft]
        cyx_tr  = cyx[ne_topright]
        cyx_bl  = cyx[ne_botleft]
        cyx_br  = cyx[ne_botright]

        cyy_tl  = cyy[ne_topleft]
        cyy_tr  = cyy[ne_topright]
        cyy_bl  = cyy[ne_botleft]
        cyy_br  = cyy[ne_botright]
        

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

        # if y_rayleigh and (cnt_y == 0):
        #     hplusx_topleft = 0.
        #     hplusy_topleft = 0. 
        #     hplusx_topright = 0.
        #     hplusy_topright = 0.

        Dx_tl = 0.5 * (topmid   - topleft + midmid   - midleft) * hplusx_topleft
        Dx_tr = 0.5 * (topright - topmid  + midright - midmid ) * hplusx_topright
        Dx_bl = 0.5 * (botmid   - botleft + midmid   - midleft) * hplusx_botleft
        Dx_br = 0.5 * (botright - botmid  + midright - midmid ) * hplusx_botright

        Dy_tl = 0.5 * (midmid   - topmid   + midleft - topleft) * hplusy_topleft
        Dy_tr = 0.5 * (midright - topright + midmid  - topmid ) * hplusy_topright
        Dy_bl = 0.5 * (botmid   - midmid   + botleft - midleft) * hplusy_botleft
        Dy_br = 0.5 * (botright - midright + botmid  - midmid ) * hplusy_botright

        fac = 1.0
        Dxx = 0.5 * (cxx_tr * Dx_tr - cxx_tl * Dx_tl + cxx_br * Dx_br - cxx_bl * Dx_bl) * oodx * oodx * fac
        Dyy = 0.5 * (cyy_br * Dy_br - cyy_tr * Dy_tr + cyy_bl * Dy_bl - cyy_tl * Dy_tl) * oody * oody * fac
        Dyx = 0.5 * (cxy_br * Dy_br - cxy_bl * Dy_bl + cxy_tr * Dy_tr - cxy_tl * Dy_tl) * oody * oodx * fac
        Dxy = 0.5 * (cyx_br * Dx_br - cyx_tr * Dx_tr + cyx_bl * Dx_bl - cyx_tl * Dx_tl) * oodx * oody * fac

        lap[idx] = Dxx + Dyy + Dyx + Dxy + hcenter[idx] * p[idx]

        lap[idx] *= diag_inv[idx]

        cnt_x += 1
        if cnt_x % iicxn == 0:
            cnt_y += 1
            cnt_x = 0
        
    return lap


@nb.stencil
def kernel_9pt(a, dx, dy, hpx, hpy, hpc, diag_inv, cxx, cyy, cxy, cyx):
    oodx = 1.0 / dx
    oody = 1.0 / dy

    topleft = a[1,-1]
    topmid = a[1,0]
    topright = a[1,1]
    
    midleft = a[0,-1]
    midmid = a[0,0]
    midright = a[0,1]
    
    botleft = a[-1,-1]
    botmid = a[-1,0]
    botright = a[-1,1]

    # shftr = 1
    # shftc = 0
    # blr = 0 - shftr
    # blc = 0 - shftc
    # brr = 0 - shftr
    # brc = 1 - shftc
    # tlr = 1 - shftr
    # tlc = 0 - shftc
    # trr = 1 - shftr
    # trc = 1 - shftc

    # hpx_bl = hpx[blr,blc]
    # hpx_br = hpx[brr,brc]
    # hpx_tl = hpx[tlr,tlc]
    # hpx_tr = hpx[trr,trc]

    # hpy_bl = hpy[blr,blc]
    # hpy_br = hpy[brr,brc]
    # hpy_tl = hpy[tlr,tlc]
    # hpy_tr = hpy[trr,trc]

    # cxx_bl = cxx[blr,blc]
    # cxx_br = cxx[brr,brc]
    # cxx_tl = cxx[tlr,tlc]
    # cxx_tr = cxx[trr,trc]

    # cyy_bl = cyy[blr,blc]
    # cyy_br = cyy[brr,brc]
    # cyy_tl = cyy[tlr,tlc]
    # cyy_tr = cyy[trr,trc]

    # cxy_bl = cxy[blr,blc]
    # cxy_br = cxy[brr,brc]
    # cxy_tl = cxy[tlr,tlc]
    # cxy_tr = cxy[trr,trc]

    # cyx_bl = cyx[blr,blc]
    # cyx_br = cyx[brr,brc]
    # cyx_tl = cyx[tlr,tlc]
    # cyx_tr = cyx[trr,trc]

    hpx_bl = hpx[0,0]
    hpx_br = hpx[0,1]
    hpx_tl = hpx[1,0]
    hpx_tr = hpx[1,1]

    hpy_bl = hpy[0,0]
    hpy_br = hpy[0,1]
    hpy_tl = hpy[1,0]
    hpy_tr = hpy[1,1]

    cxx_bl = cxx[0,0]
    cxx_br = cxx[0,1]
    cxx_tl = cxx[1,0]
    cxx_tr = cxx[1,1]

    cyy_bl = cyy[0,0]
    cyy_br = cyy[0,1]
    cyy_tl = cyy[1,0]
    cyy_tr = cyy[1,1]

    cxy_bl = cxy[0,0]
    cxy_br = cxy[0,1]
    cxy_tl = cxy[1,0]
    cxy_tr = cxy[1,1]

    cyx_bl = cyx[0,0]
    cyx_br = cyx[0,1]
    cyx_tl = cyx[1,0]
    cyx_tr = cyx[1,1]
    
    Dx_tl = 0.5 * (topmid   - topleft + midmid   - midleft) * hpx_tl
    Dx_tr = 0.5 * (topright - topmid  + midright - midmid ) * hpx_tr
    Dx_bl = 0.5 * (botmid   - botleft + midmid   - midleft) * hpx_bl
    Dx_br = 0.5 * (botright - botmid  + midright - midmid ) * hpx_br

    Dy_tl = 0.5 * (topmid   - midmid   + topleft - midleft) * hpy_tl
    Dy_tr = 0.5 * (topright - midright + topmid  - midmid ) * hpy_tr
    Dy_bl = 0.5 * (midmid   - botmid   + midleft - botleft) * hpy_bl
    Dy_br = 0.5 * (midright - botright + midmid  - botmid ) * hpy_br

    # print(hpx_tl, hpx_tr, hpx_bl, hpx_br)
    # print(hpy_tl, hpy_tr, hpy_bl, hpy_br)
    # print("")
    
    Dxx = 0.5 * (cxx_tr * Dx_tr - cxx_tl * Dx_tl + cxx_br * Dx_br - cxx_bl * Dx_bl) * oodx * oodx
    Dyy = 0.5 * (cyy_tr * Dy_tr - cyy_br * Dy_br + cyy_tl * Dy_tl - cyy_bl * Dy_bl) * oody * oody
    Dyx = 0.5 * (cxy_br * Dy_br - cxy_bl * Dy_bl + cxy_tr * Dy_tr - cxy_tl * Dy_tl) * oody * oodx
    Dxy = 0.5 * (cyx_tr * Dx_tr - cyx_br * Dx_br + cyx_tl * Dx_tl - cyx_bl * Dx_bl) * oodx * oody
    
    return ((Dxx + Dyy + Dyx + Dxy) + hpc[0,0] * a[0,0]) * diag_inv[0,0]

@nb.njit(cache=True)
def periodic(arr):
    # wall padding
    arr[0,:] = arr[2,:]
    arr[-1,:] = arr[-3,:]

    # periodic padding
    arr[:,0] = arr[:,-3]
    arr[:,-1] = arr[:,2]
    
    return arr



def stencil_27pt(elem,node,mpv,ud,diag_inv,dt):
    oodxyz = node.dxyz
    oodxyz = 1./(oodxyz**2)
    oodx2, oody2, oodz2 = oodxyz[0], oodxyz[1], oodxyz[2]
    odx, odz = 1./node.dx, 1./node.dz

    i0 = (slice(0,-1),slice(0,-1),slice(0,-1))
    i1 = (slice(1,-1),slice(1,-1),slice(1,-1))
    i2 = (slice(2,-2),slice(2,-2),slice(2,-2))

    ndim = elem.ndim
    periodicity = np.empty(ndim, dtype='int')
    for dim in range(ndim):
        periodicity[dim] = ud.bdry_type[dim] == BdryType.PERIODIC

    hplusx = mpv.wplus[0][i0][i1]
    hplusy = mpv.wplus[1][i0][i1]
    hplusz = mpv.wplus[2][i0][i1]

    hcenter = mpv.wcenter[i2]
    diag_inv = diag_inv[i1]

    corrf = dt * ud.coriolis_strength[0]

    return lambda p : lap3D(p, hplusx, hplusy, hplusz, hcenter, oodx2, oody2, oodz2, periodicity, diag_inv, corrf, odx, odz)

@nb.jit(nopython=True, cache=True, nogil=False)
def lap3D(p0, hplusx, hplusy, hplusz, hcenter, oodx2, oody2, oodz2, periodicity, diag_inv, corrf, odx, odz):
    shx, shy, shz = hcenter.shape
    p = p0.reshape(shz+2,shy+2,shx+2)

    coeff = 1./16
    lap = np.zeros_like(p)
    
    # cut out four cubes from the 3d array corresponding to the nodes... in each axial direction.
    toplefts = [(slice(0,None),slice(0,-1),slice(0,-1)),
                (slice(0,-1),slice(0,None),slice(0,-1)),
                (slice(0,-1),slice(0,-1),slice(0,None))
             ]
    toprights = [(slice(0,None),slice(0,-1),slice(1,None)),
                 (slice(1,None),slice(0,None),slice(0,-1)),
                 (slice(0,-1),slice(1,None),slice(0,None))
                ]
    
    botlefts = [(slice(0,None),slice(1,None),slice(0,-1)),
                (slice(0,-1),slice(0,None),slice(1,None)),
                (slice(1,None),slice(0,-1),slice(0,None))
               ]
    botrights = [(slice(0,None),slice(1,None),slice(1,None)),
                 (slice(1,None),slice(0,None),slice(1,None)),
                 (slice(1,None),slice(1,None),slice(0,None))
                ]

    cnt = 0
    for bc in periodicity:
        if bc == True and cnt == 0:
            tmp = p[1,:,:]
            p[0,:,:] = p[-3,:,:]
            p[-1,:,:] = p[2,:,:]
            p[1,:,:] = p[-2,:,:]
            p[-2,:,:] = tmp
        elif bc == False and cnt == 0:
            hplusx[0,:,:] = 0.0
            hplusx[-1,:,:] = 0.0
            hplusy[0,:,:] = 0.0
            hplusy[-1,:,:] = 0.0
            hplusz[0,:,:] = 0.0
            hplusz[-1,:,:] = 0.0
        if bc == True and cnt == 1:
            tmp = p[:,1,:]
            p[:,0,:] = p[:,-3,:]
            p[:,-1,:] = p[:,2,:]
            p[:,1,:] = p[:,-2,:]
            p[:,-2,:] = tmp
        elif bc == False and cnt == 1:
            hplusx[:,0,:] = 0.0
            hplusx[:,-1,:] = 0.0
            hplusy[:,0,:] = 0.0
            hplusy[:,-1,:] = 0.0
            hplusz[:,0,:] = 0.0
            hplusz[:,-1,:] = 0.0
        if bc == True and cnt == 2:
            tmp = p[:,:,1]
            p[:,:,0] = p[:,:,-3]
            p[:,:,-1] = p[:,:,2]
            p[:,:,1] = p[:,:,-2]
            p[:,:,-2] = tmp
        elif bc == False and cnt ==2:
            hplusx[:,:,0] = 0.0
            hplusx[:,:,-1] = 0.0
            hplusy[:,:,0] = 0.0
            hplusy[:,:,-1] = 0.0
            hplusz[:,:,0] = 0.0
            hplusz[:,:,-1] = 0.0
        cnt += 1

    leftz = p[:,:,:-1]
    rightz = p[:,:,1:]
    
    z_fluxes = (rightz - leftz)

    lefty = p[:,:-1,:]
    righty = p[:,1:,:]

    y_fluxes = (righty - lefty)

    leftx = p[:-1,:,:]
    rightx = p[1:,:,:]

    x_fluxes = (rightx - leftx)

    x_flx = x_fluxes[toplefts[0]] + x_fluxes[toprights[0]] + x_fluxes[botlefts[0]] + x_fluxes[botrights[0]]
    y_flx = y_fluxes[toplefts[1]] + y_fluxes[toprights[1]] + y_fluxes[botlefts[1]] + y_fluxes[botrights[1]]
    z_flx = z_fluxes[toplefts[2]] + z_fluxes[toprights[2]] + z_fluxes[botlefts[2]] + z_fluxes[botrights[2]]

    hxzp = hplusx * z_flx
    hxzpm = hxzp[:-1,:,:]
    hxzpm = hxzpm[toplefts[0]] + hxzpm[toprights[0]] + hxzpm[botlefts[0]] + hxzpm[botrights[0]]
    hxzpp = hxzp[1:,:,:]
    hxzpp = hxzpp[toplefts[0]] + hxzpp[toprights[0]] + hxzpp[botlefts[0]] + hxzpp[botrights[0]]

    hzxp = hplusz * x_flx
    hzxpm = hzxp[:,:,:-1]
    hzxpm = hzxpm[toplefts[2]] + hzxpm[toprights[2]] + hzxpm[botlefts[2]] + hzxpm[botrights[2]]
    hzxpp = hzxp[:,:,1:]
    hzxpp = hzxpp[toplefts[2]] + hzxpp[toprights[2]] + hzxpp[botlefts[2]] + hzxpp[botrights[2]]

    x_flx = hplusx * x_flx
    x_flxm = x_flx[:-1,:,:]
    x_flxm = x_flxm[toplefts[0]] + x_flxm[toprights[0]] + x_flxm[botlefts[0]] + x_flxm[botrights[0]]
    x_flxp = x_flx[1:,:,:]
    x_flxp = x_flxp[toplefts[0]] + x_flxp[toprights[0]] + x_flxp[botlefts[0]] + x_flxp[botrights[0]]

    y_flx = hplusy * y_flx
    y_flxm = y_flx[:,:-1,:]
    y_flxm = y_flxm[toplefts[1]] + y_flxm[toprights[1]] + y_flxm[botlefts[1]] + y_flxm[botrights[1]]
    y_flxp = y_flx[:,1:,:]
    y_flxp = y_flxp[toplefts[1]] + y_flxp[toprights[1]] + y_flxp[botlefts[1]] + y_flxp[botrights[1]]

    z_flx = hplusz * z_flx
    z_flxm = z_flx[:,:,:-1]
    z_flxm = z_flxm[toplefts[2]] + z_flxm[toprights[2]] + z_flxm[botlefts[2]] + z_flxm[botrights[2]]
    z_flxp = z_flx[:,:,1:]
    z_flxp = z_flxp[toplefts[2]] + z_flxp[toprights[2]] + z_flxp[botlefts[2]] + z_flxp[botrights[2]]

    lap[1:-1,1:-1,1:-1] = oodx2 * coeff * (-x_flxm + x_flxp) + \
                          oody2 * coeff * (-y_flxm + y_flxp) + \
                          oodz2 * coeff * (-z_flxm + z_flxp) + \
                          +1.0 * odx * odz * coeff * corrf * (hxzpp - hxzpm) + \
                          -1.0 * odx * odz * coeff * corrf * (hzxpp - hzxpm) + \
                          hcenter * p[1:-1,1:-1,1:-1]

    lap = lap * diag_inv

    return lap



def stencil_hs(elem,node,mpv,ud,diag_inv,dt):
    oodx2 = 1./node.dx**2
    oodz2 = 1./node.dz**2
    odx, odz = 1./node.dx, 1./node.dz
    igy = elem.igs[1]

    proj = (slice(None,),igy,slice(None,))
    i0 = (slice(0,-1),slice(0,-1))
    i1 = (slice(1,-1),slice(1,-1))
    i2 = (slice(2,-2),slice(2,-2))

    ndim = elem.ndim
    periodicity = np.empty(ndim, dtype='int')
    for dim in range(0,ndim,2):
        periodicity[dim] = ud.bdry_type[dim] == BdryType.PERIODIC

    hplusx = mpv.wplus[0][proj][i0][i1]
    hplusz = mpv.wplus[2][proj][i0][i1]
    hcenter = mpv.wcenter[proj][i2]
    diag_inv = diag_inv[proj][i1]

    corrf = dt * ud.coriolis_strength[1]

    return lambda p : lapHS(p, hplusx, hplusz, hcenter, oodx2, oodz2, periodicity, diag_inv, corrf, odx, odz)


@nb.jit(nopython=True, cache=True, nogil=True)
def lapHS(p0, hplusx, hplusz, hcenter, oodx2, oodz2, periodicity, diag_inv, corrf, odx, odz):
    shx, shz = hcenter.shape
    p = p0.reshape(shx+2,shz+2)

    coeff = 1./4
    lap = np.zeros_like(p)
    
    cnt = 0
    for bc in periodicity:
        if bc == True and cnt == 0:
            tmp = p[1,:]
            p[0,:] = p[-3,:]
            p[-1,:] = p[2,:]
            p[1,:] = p[-2,:]
            p[-2,:] = tmp

            # p[0,:] = p[-4,:]
            # p[-1,:] = p[3,:]
            # p[1,:] = p[-3,:]
            # p[-2,:] = p[2,:]
        elif bc == False and cnt == 0:
            hplusx[0,:] = 0.0
            hplusx[-1,:] = 0.0
            hplusz[0,:] = 0.0
            hplusz[-1,:] = 0.0
        if bc == True and cnt == 2:
            tmp = p[:,1]
            p[:,0] = p[:,-3]
            p[:,-1] = p[:,2]
            p[:,1] = p[:,-2]
            p[:,-2] = tmp
            # p[:,0] = p[:,-4]
            # p[:,-1] = p[:,3]
            # p[:,1] = p[:,-3]
            # p[:,-2] = p[:,2]
        elif bc == False and cnt == 2:
            hplusx[:,0] = 0.0
            hplusx[:,-1] = 0.0
            hplusz[:,0] = 0.0
            hplusz[:,-1] = 0.0
        cnt += 1
    
    leftz = p[:,:-1]
    rightz = p[:,1:]
    
    z_fluxes = (rightz - leftz)

    leftx = p[:-1,:]
    rightx = p[1:,:]

    x_fluxes = (rightx - leftx)

    x_flx = x_fluxes[:,:-1] + x_fluxes[:,1:]
    z_flx = z_fluxes[:-1,:] + z_fluxes[1:,:]

    hxzp = hplusx * z_flx
    hxzpm = hxzp[:-1,:]
    hxzpm = hxzpm[:,:-1] + hxzpm[:,1:]
    # hxzpm = hxzpm[:-1,:] + hxzpm[1:,:]
    hxzpp = hxzp[1:,:]
    hxzpp = hxzpp[:,:-1] + hxzpp[:,1:]
    # hxzpp = hxzpp[:-1,:] + hxzpp[1:,:]

    hzxp = hplusz * x_flx
    hzxpm = hzxp[:,:-1]
    hzxpm = hzxpm[:-1,:] + hzxpm[1:,:]
    # hzxpm = hzxpm[:,:-1] + hzxpm[:,1:]
    hzxpp = hzxp[:,1:]
    hzxpp = hzxpp[:-1,:] + hzxpp[1:,:]
    # hzxpp = hzxpp[:,:-1] + hzxpp[:,1:]

    x_flx = hplusx * x_flx
    x_flxm = x_flx[:-1,:]
    x_flxm = x_flxm[:,:-1] + x_flxm[:,1:]
    x_flxp = x_flx[1:,:]
    x_flxp = x_flxp[:,:-1] + x_flxp[:,1:]

    z_flx = hplusz * z_flx
    z_flxm = z_flx[:,:-1]
    z_flxm = z_flxm[:-1,:] + z_flxm[1:,:]
    z_flxp = z_flx[:,1:]
    z_flxp = z_flxp[:-1,:] + z_flxp[1:,:]

    lap[1:-1,1:-1] = oodx2 * coeff * (-x_flxm + x_flxp) + \
                     oodz2 * coeff * (-z_flxm + z_flxp) + \
                     +1.0 * odx * odz * coeff * corrf * (hxzpp - hxzpm) + \
                     -1.0 * odx * odz * coeff * corrf * (hzxpp - hzxpm) + \
                     hcenter * p[1:-1,1:-1]

    lap = lap * diag_inv

    return lap



def stencil_vs(elem,node,mpv,ud,diag_inv,dt):
    igx = elem.igx
    igy = elem.igy
    igz = elem.igz

    icxn = node.icx
    icyn = node.icy

    iicxn = icxn - (2 * igx)
    iicyn = icyn - (2 * igy)

    iicxn, iicyn = iicyn, iicxn

    dx = node.dy
    dy = node.dx

    oodx2 = 0.5 / (dx**2)
    oody2 = 0.5 / (dy**2)

    xx, yy = 1, 0

    proj = (slice(None,),slice(None,),igz)
    i0 = (slice(0,-1),slice(0,-1))
    i1 = (slice(1,-1),slice(1,-1))
    i2 = (slice(igx,-igx),slice(igy,-igy))

    VS = True

    if not VS:
        x_periodic = ud.bdry_type[xx] == BdryType.PERIODIC
        y_periodic = ud.bdry_type[yy] == BdryType.PERIODIC

        x_wall = ud.bdry_type[xx] == BdryType.WALL
        y_wall = ud.bdry_type[yy] == BdryType.WALL

        hplusx = mpv.wplus[xx][proj][i2].reshape(-1,)
        hplusy = mpv.wplus[yy][proj][i2].reshape(-1,)
        hcenter = mpv.wcenter[proj][i2].reshape(-1,)
        diag_inv = diag_inv[proj][i2].reshape(-1,)

        return lambda p : lap2D(p, igx,igy, iicxn, iicyn, hplusx, hplusy, hcenter, oodx2, oody2, x_periodic, y_periodic, x_wall, y_wall, diag_inv)

    if VS:
        oodx2 = 1./node.dx**2
        oody2 = 1./node.dy**2

        ndim = elem.ndim
        periodicity = np.zeros(ndim, dtype='int')
        for dim in range(0,ndim,2):
            periodicity[dim] = ud.bdry_type[dim] == BdryType.PERIODIC

        hplusx = mpv.wplus[0][proj][i0][i1]
        hplusy = mpv.wplus[1][proj][i0][i1]
        hcenter = mpv.wcenter[proj][i2]
        diag_inv = diag_inv[proj][i1]

        return lambda p : lapVS(p, hplusx, hplusy, hcenter, oodx2, oody2, periodicity, diag_inv)


@nb.jit(nopython=True, cache=True, nogil=False)
def lapVS(p0, hplusx, hplusy, hcenter, oodx2, oody2, periodicity, diag_inv):
    shx, shy = hcenter.shape
    p = p0.reshape(shx+2,shy+2)

    coeff = 1./4
    lap = np.zeros_like(p)
    
    cnt = 0
    for bc in periodicity:
        if bc == True and cnt == 0:
            p[0,:] = p[-3,:]
            p[-1,:] = p[2,:]

        elif bc == False and cnt == 0:
            hplusx[0,:] = 0.0
            hplusx[-1,:] = 0.0
            hplusy[0,:] = 0.0
            hplusy[-1,:] = 0.0

            # p[0,:] = p[2,:]
            # p[-1,:] = p[-3,:]

        if bc == True and cnt == 1:
            p[:,0] = p[:,-3]
            p[:,-1] = p[:,2]

        elif bc == False and cnt == 1:
            hplusx[:,0] = 0.0
            hplusx[:,-1] = 0.0
            hplusy[:,0] = 0.0
            hplusy[:,-1] = 0.0

            # p[:,0] = p[:,2]
            # p[:,-1] = p[:,-3]
            
        cnt += 1
    
    lefty = p[:,:-1]
    righty = p[:,1:]
    
    y_fluxes = (righty - lefty)

    leftx = p[:-1,:]
    rightx = p[1:,:]

    x_fluxes = (rightx - leftx)

    x_flx = x_fluxes[:,:-1] + x_fluxes[:,1:]
    y_flx = y_fluxes[:-1,:] + y_fluxes[1:,:]

    x_flx = hplusx * x_flx
    x_flxm = x_flx[:-1,:]
    x_flxm = x_flxm[:,:-1] + x_flxm[:,1:]
    x_flxp = x_flx[1:,:]
    x_flxp = x_flxp[:,:-1] + x_flxp[:,1:]

    y_flx = hplusy * y_flx
    y_flxm = y_flx[:,:-1]
    y_flxm = y_flxm[:-1,:] + y_flxm[1:,:]
    y_flxp = y_flx[:,1:]
    y_flxp = y_flxp[:-1,:] + y_flxp[1:,:]

    lap[1:-1,1:-1] = oodx2 * coeff * (-x_flxm + x_flxp) + \
                     oody2 * coeff * (-y_flxm + y_flxp) + \
                     hcenter * p[1:-1,1:-1]
                     
    lap = lap * diag_inv

    return lap



def precon_diag_prepare(mpv, elem, node, ud, coriolis):
    dx, dy, dz = node.dx, node.dy, node.dz

    x_periodic = ud.bdry_type[0] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    z_periodic = ud.bdry_type[2] == BdryType.PERIODIC
    periodicity = (x_periodic, y_periodic, z_periodic)
    # igx = node.igx
    # igy = node.igy
    igs = node.igs
    ndim = elem.ndim

    # idx_e = (slice(igx - x_periodic,-igx + x_periodic - 1),slice(igy - y_periodic, -igy + y_periodic - 1))
    idx_periodic = [slice(None)] * elem.ndim
    idx_n, idx_e = np.copy(idx_periodic), np.copy(idx_periodic)

    for dim in range(ndim):
        if ud.bdry_type[dim] == BdryType.PERIODIC:
            idx_periodic[dim] = slice(1,-1)

        idx_e[dim] = slice(igs[dim] - periodicity[dim], -igs[dim] + periodicity[dim] - 1)
        idx_n[dim] = slice(igs[dim], -igs[dim])

    idx_periodic, idx_e, idx_n = tuple(idx_periodic), tuple(idx_e), tuple(idx_n)

    # idx_n = (slice(igx,-igx), slice(igy,-igy))

    hplusxx = mpv.wplus[0] #* (coriolis[0].T)
    hplusyy = mpv.wplus[1] #* (coriolis[1].T)
    hplusxy = mpv.wplus[0] #* (coriolis[3].T)
    hplusyx = mpv.wplus[1] #* (coriolis[2].T)

    if ndim == 3: hplusz = mpv.wplus[2]

    nine_pt = 0.25 * (2.0) * 1.0
    nine_pt = 0.5 * 0.5 
    if ndim == 2:
        coeff = 1.0 - nine_pt
    elif ndim == 3:
        coeff = 0.0625
    else:
        assert 0, "ndim = 1?"

    wxx = coeff / (dx**2)
    wyy = coeff / (dy**2)
    wzz = coeff / (dz**2)

    wxy = coeff / (dx * dy)
    wyx = coeff / (dy * dx)

    diag_kernel = np.array(np.ones([2] * ndim))

    # diag = np.zeros((node.sc)).squeeze()
    # diag[idx_n] = -wx * signal.fftconvolve(hplusx[idx_e],diag_kernel,mode='full')[idx_periodic] 
    # diag[idx_n] -= wy * signal.fftconvolve(hplusy[idx_e],diag_kernel,mode='full')[idx_periodic]
    # if ndim == 3:
    #     diag[idx_n] -= wz * signal.fftconvolve(hplusz[idx_e],diag_kernel,mode='full')[idx_periodic]

    # diag[idx_n] += mpv.wcenter[idx_n]
    # diag[idx_n] = 1.0 / diag[idx_n]

    diag = np.zeros_like(mpv.wcenter)
    diag[...] = -wxx * signal.fftconvolve(hplusxx,diag_kernel,mode='valid')
    diag[...] -= wyy * signal.fftconvolve(hplusyy,diag_kernel,mode='valid')
    diag[...] -= wxy * signal.fftconvolve(hplusxy,diag_kernel,mode='valid')
    diag[...] -= wyx * signal.fftconvolve(hplusyx,diag_kernel,mode='valid')
    if ndim == 3:
        diag[...] -= wzz * signal.fftconvolve(hplusz,diag_kernel,mode='valid')

    diag[...] += mpv.wcenter
    diag[...] = 1.0 / diag

    return diag





def stencil_3pt(elem,node,ud):
    # 1d stencil for exner pressure perturbation constraint in 2d
    igx = elem.igx
    igy = elem.igy

    icxn = node.icx
    icyn = node.icy

    iicxn = icxn - (2 * igx)
    iicyn = icyn - (2 * igy)

    iicxn, iicyn = iicyn, iicxn

    dx = node.dx
    # print(node.dx)

    oodx2 = 1. / (dx**2)

    # x_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    x_periodic = ud.bdry_type[0] == BdryType.PERIODIC
    x_wall = ud.bdry_type[0] == BdryType.WALL

    return lambda p : lap2D_exner(p,iicxn, iicyn, oodx2, x_periodic, x_wall)

@jit(nopython=True, nogil=True, cache=True)
def lap2D_exner(p, iicxn, iicyn, oodx2, x_periodic, x_wall):
    ngnc = (iicxn) * (iicyn)
    lap = np.zeros((ngnc))
    cnt_x = 0
    cnt_y = 0

    for idx in range(iicxn * iicyn):
        # left_idx = idx - 1
        # mid_idx = idx
        # right_idx = idx + 1

        # if cnt_x == 0:
        #     left_idx += iicxn - 1

        #     if x_periodic:
        #         mid_idx += iicxn - 1

        # if cnt_x == (iicxn - 1):
        #     right_idx -= iicxn - 1

        #     if x_periodic:
        #         mid_idx -= iicxn - 1

        right_idx = idx - iicxn
        mid_idx = idx
        left_idx = idx + iicxn

        if cnt_y == 0:
            if x_periodic:
                right_idx += ((iicxn) * (iicyn - 1))
                mid_idx += ((iicxn) * (iicyn - 1))

        if cnt_y == (iicyn - 1):
            if x_periodic:
                left_idx -= ((iicxn) * (iicyn - 1))
                mid_idx -= ((iicxn) * (iicyn - 1))

        left = p[left_idx]
        mid = p[mid_idx]
        right = p[right_idx]

        lap[idx] = oodx2 * (left - 2. * mid + right)

        cnt_x += 1
        if cnt_x % iicxn == 0:
            cnt_y += 1
            cnt_x = 0
        
    return lap


def stencil_5pt(elem,node,ud):
    igx = elem.igx
    igy = elem.igy

    icxn = node.icx
    icyn = node.icy

    iicxn = icxn - (2 * igx)
    iicyn = icyn - (2 * igy)

    iicxn, iicyn = iicyn, iicxn

    dx = node.dy
    dy = node.dx

    oodx2 = 1.0 / (dx**2)
    oody2 = 1.0 / (dy**2)

    # oodx2 = 1.0
    # oody2 = 1.0

    x_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[0] == BdryType.PERIODIC

    x_wall = ud.bdry_type[1] == BdryType.WALL
    y_wall = ud.bdry_type[0] == BdryType.WALL

    return lambda p : lap2D_5pt(p, iicxn, iicyn, oodx2, oody2, x_periodic, y_periodic, x_wall, y_wall)

@jit(nopython=True, nogil=True, cache=True)
def lap2D_5pt(p, iicxn, iicyn, oodx2, oody2, x_periodic, y_periodic, x_wall, y_wall):
    ngnc = (iicxn) * (iicyn)
    lap = np.zeros((ngnc))
    cnt_x = 0
    cnt_y = 0

    for idx in range(iicxn * iicyn):

        # get indices of the 9pt stencil
        midleft_idx = idx - 1
        topmid_idx = idx - iicxn
        midmid_idx = idx
        botmid_idx = idx + iicxn
        midright_idx = idx + 1

        if cnt_x == 0:
            midleft_idx += iicxn - 1
            if x_periodic:
                topmid_idx += iicxn - 1
                midmid_idx += iicxn - 1
                botmid_idx += iicxn - 1

        if cnt_x == (iicxn - 1):
            midright_idx -= iicxn - 1

            if x_periodic:
                topmid_idx -= iicxn - 1
                midmid_idx -= iicxn - 1
                botmid_idx -= iicxn - 1

        if cnt_y == 0:
            topmid_idx += ((iicxn) * (iicyn - 1))

            if y_periodic:
                midleft_idx += ((iicxn) * (iicyn - 1))
                midmid_idx += ((iicxn) * (iicyn - 1))
                midright_idx += ((iicxn) * (iicyn - 1))

        if cnt_y == (iicyn - 1):
            botmid_idx -= ((iicxn) * (iicyn - 1))

            if y_periodic:
                midleft_idx -= ((iicxn) * (iicyn - 1))
                midmid_idx -= ((iicxn) * (iicyn - 1))
                midright_idx -= ((iicxn) * (iicyn - 1))

        midleft = p[midleft_idx]
        topmid = p[topmid_idx]
        midmid = p[midmid_idx]
        botmid = p[botmid_idx]
        midright = p[midright_idx]

        if x_wall and (cnt_x == 0):
            midleft = 0.

        if x_wall and (cnt_x == (iicxn - 1)):
            midright = 0.

        if y_wall and (cnt_y == 0):
            topmid = 0.
            
        if y_wall and (cnt_y == (iicyn - 1)):
            botmid = 0.
                                
        lap[idx] = oodx2 * (midleft - 2.0 * midmid + midright) + oody2 * (topmid - 2.0 * midmid + botmid)

        cnt_x += 1
        if cnt_x % iicxn == 0:
            cnt_y += 1
            cnt_x = 0
        
    return lap