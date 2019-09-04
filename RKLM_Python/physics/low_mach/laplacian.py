import numpy as np
from inputs.enum_bdry import BdryType
from scipy import sparse, signal
from itertools import product
from numba import jit

def stencil_9pt_operator(elem,node,mpv,ud):
    ndim = node.ndim
    igs = node.igs

    inner_idx = np.empty((ndim), dtype=object)
    inner_nidx_periodic = np.empty_like(inner_idx)
    inner_eidx_periodic = np.empty_like(inner_idx)
    for dim in range(ndim):
        is_periodic = ud.bdry_type[dim] == BdryType.PERIODIC
        inner_nidx_periodic[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic)
        inner_eidx_periodic[dim] = slice(igs[dim]-is_periodic,-igs[dim]+is_periodic-1)
        inner_idx[dim] = slice(igs[dim],-igs[dim])

    inner_idx, inner_nidx_periodic = tuple(inner_idx), tuple(inner_nidx_periodic)
    inner_eidx_periodic = tuple(inner_eidx_periodic)
    four_quads_indices = [idx for idx in product([slice(0,-1),slice(1,None)], repeat=ndim)]

    hplusx, hplusy = mpv.wplus[0], mpv.wplus[1]

    xflux_mid_mid = np.zeros_like(hplusx[inner_idx])
    for quad in four_quads_indices:
        xflux_mid_mid += hplusx[inner_eidx_periodic][quad]
    
    xflux_mid_mid *= -0.5
    xflux_mid_left = hplusx[inner_eidx_periodic][:-1,-1] + hplusx[inner_eidx_periodic][:-1,1:] 
    xflux_mid_left *= 0.5
    xflux_mid_right = hplusx[inner_eidx_periodic][1:,-1] + hplusx[inner_eidx_periodic][1:,1:] 
    xflux_mid_right *= 0.5

    tmp = np.diag(xflux_mid_mid[1:,1:].ravel())
    diagonals_mid = [xflux_mid_left.ravel()[:-1], xflux_mid_mid.ravel(), xflux_mid_right.ravel()]
    offsets = [-1,0,1]
    tmp = sparse.diags(diagonals_mid,offsets,format='lil')

    # print(xflux_mid_right.shape)



def stencil_9pt(elem,node,mpv,ud):
    igx = elem.igx
    igy = elem.igy

    icx = elem.icx
    icxn = node.icx
    icy = elem.icy
    icyn = node.icy

    nc = node.sc

    dx = node.dy
    dy = node.dy

    hplusx = mpv.wplus[0].reshape(-1,)
    hplusy = mpv.wplus[1].reshape(-1,)
    hcenter = mpv.wcenter.reshape(-1,)

    oodx2 = 0.5 / (dx**2)
    oody2 = 0.5 / (dy**2)
    nine_pt = 0.25 * (2.0) * 1.0

    x_periodic = ud.bdry_type[0] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[1] == BdryType.PERIODIC

    lap = np.zeros((53**2))
    def lap2D(p):
        # lap = np.zeros_like(p)

        for j in range(igy - y_periodic, icy -igy + y_periodic):
            me = j * icx
            mn = j * icxn

            for i in range(igx - x_periodic, icx -igx + x_periodic):
                ne = me + i
                nn = mn + i
                nn1 = nn + 1
                nnicxn = nn + icxn
                nn1icxn = nn + 1 + icxn

                dsq_p_dxdy = p[nn1icxn] - p[nnicxn] - p[nn1] + p[nn]

                flux_x_lower  = hplusx[ne] * oodx2 * ( (p[nn1]     - p[nn]    ) + nine_pt * dsq_p_dxdy)
                flux_x_upper  = hplusx[ne] * oodx2 * ( (p[nn1icxn] - p[nnicxn]) - nine_pt * dsq_p_dxdy)
                
                flux_y_left   = hplusy[ne] * oody2 * ( (p[nnicxn]  - p[nn]    ) + nine_pt * dsq_p_dxdy)
                flux_y_right  = hplusy[ne] * oody2 * ( (p[nn1icxn] - p[nn1]   ) - nine_pt * dsq_p_dxdy)
                
                lap[nn]      += (  flux_x_lower + flux_y_left )
                lap[nn1]     += (- flux_x_lower + flux_y_right)
                lap[nnicxn]  += (  flux_x_upper - flux_y_left )
                lap[nn1icxn] += (- flux_x_upper - flux_y_right)

        for j in range(igy,icyn-igy):
            mn = j * icxn
            for i in range(igx,icxn-igx):
                nn = mn + i
                lap[nn] += hcenter[nn] * p[nn]
        return lap

    return lap2D
    

def stencil_9pt_2nd_try(rhs,elem,node,mpv,ud):
    igx = elem.igx
    igy = elem.igy

    icx = elem.icx
    icxn = node.icx
    icy = elem.icy
    icyn = node.icy

    sc = node.sc

    iicxn = icxn - (2 * igx)
    # iicxn = icxn - (1 * igx)
    iicx = icx - (2 * igx)
    iicyn = icyn - (2 * igy)
    # iicyn = icyn - (1 * igy)
    iicy = icy - (2 *igy)
    ngnc = (iicxn) * (iicyn)

    # iicxn, iicyn = iicyn, iicxn
    dx = node.dx
    dy = node.dy

    inner_domain = (slice(igx,-igx),slice(igy,-igy))
    # hplusx = mpv.wplus[0][inner_domain].reshape(-1,)
    # hplusy = mpv.wplus[1][inner_domain].reshape(-1,)
    hplusx = mpv.wplus[0][igx:-igx,igy:-igy].reshape(-1,)
    hplusy = mpv.wplus[1][igx:-igx,igy:-igy].reshape(-1,)
    # hplusx = mpv.wplus[0][igx-1:-igx+1,igy-1:-igy+1].reshape(-1,)
    # hplusy = mpv.wplus[1][igx-1:-igx+1,igy-1:-igy+1].reshape(-1,)
    # hplusx = np.vstack((hplusx,hplusx[0],hplusx[1])).reshape(-1,)

    hcenter = mpv.wcenter[inner_domain].reshape(-1,)
    # hcenter = mpv.wcenter.reshape(-1,)
    # hplusx = mpv.wplus[0].reshape(-1,)
    # hplusy = mpv.wplus[1].reshape(-1,)

    oodx2 = 0.5 / (dx**2)
    oody2 = 0.5 / (dy**2)
    nine_pt = 0.25 * (2.0) * 1.0

    x_periodic = ud.bdry_type[0] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    x_wall = ud.bdry_type[0] == BdryType.WALL
    y_wall = ud.bdry_type[1] == BdryType.WALL
    # print(x_wall)
    # print(y_wall)
    lap = np.zeros((ngnc))
    # print(lap.shape)
    # print(hcenter.shape)

    def lap2D_try(p):
        # lap = np.zeros_like(p)
        cnt_x = 1
        cnt_y = 0
        idx = 0
        cnt = 0
        cnt_idx = 0

        for idx in range(iicxn * iicyn):
            # update counter for current row, column
            if cnt_x % iicxn == 1:
                cnt_x = 1
                cnt_y += 1
            
            ne_topleft = idx - iicx - 1
            ne_topright = idx - iicx
            ne_botleft = idx - 1
            ne_botright = idx

            # get indices of the 9pt stencil
            topleft = idx - iicxn - 1
            midleft = idx - 1
            botleft = idx + iicxn - 1

            topmid = idx - iicxn
            midmid = idx
            botmid = idx + iicxn

            topright = idx - iicxn + 1
            midright = idx + 1
            botright = idx + iicxn + 1

            # do periodic boundary checks
            # x-axis periodic:
            # if x_periodic:
            topx = 1
            botx = iicxn 
            topy = 1
            boty = iicyn 

            if cnt_x == topx:
                # topleft = idx - 1
                # midleft = idx + iicxn - 1
                # botleft = idx + 2 * iicxn - 1
                
                # shift = -1
                # topmid += (iicxn + shift)
                # midmid += (iicxn + shift)
                # botmid += (iicxn + shift)

                shift = -1
                topleft += (iicxn + shift)
                midleft += (iicxn + shift)
                botleft += (iicxn + shift)

                # if 50 < idx < 100:
                #     print("cnt_x == 1, topleft =", cnt_x, cnt_y, topleft, idx)
                #     print("cnt_x == 1, midleft =", cnt_x, cnt_y, midleft, idx)
                #     print("cnt_x == 1, botleft =", cnt_x, cnt_y, botleft, idx)
            
            if cnt_x == botx:
                # topright = idx - 2 * iicxn + 1
                # midright = idx - icxn +1
                # botright = idx + 1

                # shift = -1
                # topmid -= (iicxn + shift)
                # midmid -= (iicxn + shift)
                # botmid -= (iicxn + shift)
                
                shift = -1
                topright -= (iicxn + shift)
                midright -= (iicxn + shift)
                botright -= (iicxn + shift)
                
                # if  50 < idx < 100:
                #     print("cnt_x == 1, topright =", cnt_x, cnt_y, topright, idx)
                #     print("cnt_x == 1, midright =", cnt_x, cnt_y, midright, idx)
                #     print("cnt_x == 1, botright =", cnt_x, cnt_y, botright, idx)

            # if cnt_x < 100 and idx < 100:
            #     print(cnt_x, idx)
            if x_periodic:
                if cnt_x == topx:
                    shift = 1
                    ne_topleft += (iicxn - shift)
                    ne_botleft += (iicxn - shift)

                if cnt_x == botx:
                    shift = 1
                    ne_topright -= (iicxn - shift)
                    ne_botright -= (iicxn - shift)
        
                # if cnt_x == (iicxn - 1) and (50 < ne < 100):
                #     print(ne, ne_topleft, ne_topright, ne_botleft, ne_botright)
    
                # if cnt_x == (iicxn) and (ne < 50):
                #     print(ne, ne_topleft, ne_topright, ne_botleft, ne_botright)
    
            if cnt_y == topy:
                # topleft = idx + (iicxn * (iicyn - 1)) - 1
                # topmid = idx + iicxn * (iicyn - 1)
                # topright = idx + iicxn * (iicyn - 1) + 1

                # shift = 0
                # midleft += (iicxn * (iicyn - shift))
                # midmid += (iicxn * (iicyn - shift))
                # midright += (iicxn * (iicyn - shift))

                shift = 1
                topleft += (iicxn * (iicyn - shift))
                topmid += (iicxn * (iicyn - shift))
                topright += (iicxn * (iicyn - shift))
                
                # if  0 < idx < 50:
                #     print("cnt_y == 1, topleft =", cnt_x, cnt_y, topleft, idx)
                #     print("cnt_y == 1, topmid =", cnt_x, cnt_y, topmid, idx)
                #     print("cnt_y == 1, topright =", cnt_x, cnt_y, topright, idx)

            if cnt_y == boty:
                # print("botleft", botleft, idx, iicxn, iicyn, idx - iicxn * (iicyn - 1) - 1)
                # botleft = idx - iicxn * (iicyn - 1) - 1
                # botmid = idx - iicxn * (iicyn - 1)
                # botright = idx - iicxn * (iicyn - 1) + 1

                # shift = 1
                # midleft -= (iicxn * (iicyn - shift))
                # midmid -= (iicxn * (iicyn - shift))
                # midright -= (iicxn * (iicyn - shift))

                shift = 1
                botleft -= (iicxn * (iicyn - shift))
                botmid -= (iicxn * (iicyn - shift))
                botright -= (iicxn * (iicyn - shift))

                
                # print("cnt_y == 1, topleft =", cnt_x, cnt_y, topleft, idx)
                # print("cnt_y == 1, topmid =", cnt_x, cnt_y, topmid, idx)
                # print("cnt_y == 1, topright =", cnt_x, cnt_y, topright, idx)
            if y_periodic:
                if cnt_y == topy:
                    shift = 1
                    ne_topright += ((iicyn - shift) * iicxn)
                    ne_topleft += ((iicyn - shift) * iicxn)

                if cnt_y == boty:
                    shift = 1
                    ne_botleft -= ((iicyn - shift) * iicxn)
                    ne_botright -= ((iicyn - shift) * iicxn)

            # if ne_botleft == -1:
            #     ne_botleft += 1
            # if ne_botright == 0:
            #     ne_botright += 1

            # if cnt_y == (iicyn):
            #     print(ne, ne_topleft, ne_topright, ne_botleft, ne_botright)

            # if cnt_x == (iicxn - 1) and idx < 100:
            #     print("iicxn - 1: ", ne, ne_topleft, ne_topright, ne_botleft, ne_botright)

            # if cnt_x == iicxn and idx < 100:
            #     print("iicxn: ", ne, ne_topleft, ne_topright, ne_botleft, ne_botright)

            # if cnt_x == 1 and idx < 100:
            #     print("1: ", ne, ne_topleft, ne_topright, ne_botleft, ne_botright)

            # do wall-boundary checks
            # if x_wall:
            #     if cnt_x == 1:
            #         shift = 1
            #         topleft += (shift)
            #         midleft += (shift)
            #         botleft += (shift)

            #     if cnt_x == iicxn:
            #         shift = 1
            #         topright -= (shift)
            #         midright -= (shift)
            #         botright -= (shift)
                
            #     if cnt_x == 1:
            #         shift = 1
            #         ne_topleft += (shift)
            #         ne_botleft += (shift)

            #     if cnt_x == iicxn:
            #         shift = 1
            #         ne_topright -= (shift)
            #         ne_botright -= (shift)

            # if y_wall:
            #     if cnt_y == 1:
            #         shift = 1
            #         topleft += (iicxn * (shift))
            #         topmid += (iicxn * (shift))
            #         topright += (iicxn * (shift))

            #     if cnt_y == iicyn:
            #         shift = 1
            #         botleft -= (iicxn * (shift))
            #         botmid -= (iicxn * (shift))
            #         botright -= (iicxn * (shift))

            #     if cnt_y == 1:
            #         shift = 1
            #         ne_topright += ((shift) * iicxn)
            #         ne_topleft += ((shift) * iicxn)

            #     if cnt_y == iicyn:
            #         shift = 1
            #         ne_botleft -= ((shift) * iicxn)
            #         ne_botright -= ((shift) * iicxn)
            
            # get values at indices of the 9pt stencil
            nine_pt = 0.25 * (2.0) * 1.0

            oodx2 = 0.5 / (dx**2)
            oody2 = 0.5 / (dy**2)

            topleft = p[topleft]
            midleft = p[midleft]
            botleft = p[botleft]

            topmid = p[topmid]
            midmid = p[midmid]
            botmid = p[botmid]

            topright = p[topright]
            midright = p[midright]
            botright = p[botright]

            hplusx_topleft = hplusx[ne_topleft]
            hplusy_topleft = hplusy[ne_topleft]
            hplusx_topright = hplusx[ne_topright]
            hplusy_topright = hplusy[ne_topright]
            hplusx_botleft = hplusx[ne_botleft]
            hplusy_botleft = hplusy[ne_botleft]
            hplusx_botright = hplusx[ne_botright]
            hplusy_botright = hplusy[ne_botright]

            hcenter_idx = hcenter[idx]
            
            c1inx = midmid - midleft
            c1iny = midmid - topmid

            c2inx = midright - midmid
            c2iny = midmid - topmid

            c3inx = midmid - midleft
            c3iny = botmid - midmid

            c4inx = midright - midmid
            c4iny = botmid - midmid

            # let's check wall bcs
            if x_wall:
                if cnt_x == 1:
                    hplusx_topleft = 0.
                    hplusy_topleft = 0.
                    hplusx_botleft = 0.
                    hplusy_botleft = 0.
                    
                    # topleft = 0.
                    # midleft = 0.
                    # botleft = 0.
                if cnt_x == iicxn:
                # if cnt_y == iicyn:
                    hplusx_topright = 0.
                    hplusy_topright = 0.
                    hplusx_botright = 0.
                    hplusy_botright = 0.

                    # topright = 0.
                    # midright = 0.
                    # botright = 0.

            if y_wall:
                coeff = 0.
                # print(True)
                # c1inx = 0.
                # c2inx = 0.
                # c3inx = 0.
                # c4inx = 0.

                if cnt_y == topy:
                    cnt += 1
                    # None
            # if x_wall:
            #     if cnt_x == 1:
                    # print("y_wall true")
                    hplusx_topleft = coeff * hplusx_botleft
                    hplusy_topleft = coeff * hplusy_botleft
                    hplusx_topright = coeff * hplusx_botright
                    hplusy_topright = coeff * hplusy_botright

                    # hcenter_idx *= 2.
                    # nine_pt *= 2.
                    # topleft = coeff*midleft
                    # topmid = coeff*midmid
                    # topright = coeff*midright
                    # oodx2 *= 0.5
                    # oody2 *= 0.5

                if cnt_y == boty:
                    cnt += 1
                    # None
                    hplusx_botleft = coeff * hplusx_topleft
                    hplusy_botleft = coeff * hplusy_topleft
                    hplusx_botright = coeff * hplusx_topright
                    hplusy_botright = coeff * hplusy_topright

                    # hcenter_idx *= 2.
                    # nine_pt *= 2.
                    # oodx2 *= 0.5
                    # oody2 *= 0.5


            dp2dxdy1 = midmid - midleft - topmid + topleft
            dp2dxdy1 *= nine_pt
            dp2dxdy2 = midright - midmid - topright + topmid
            dp2dxdy2 *= nine_pt
            dp2dxdy3 = botmid - botleft - midmid + midleft
            dp2dxdy3 *= nine_pt
            dp2dxdy4 = botright - botmid - midright + midmid
            dp2dxdy4 *= nine_pt

            lap[idx] = - hplusx_topleft * oodx2 * (c1inx - dp2dxdy1) \
                    -  hplusy_topleft * oody2 * (c1iny - dp2dxdy1) \
                    +  hplusx_topright * oodx2 * (c2inx - dp2dxdy2) \
                    -  hplusy_topright * oody2 * (c2iny + dp2dxdy2) \
                    -  hplusx_botleft * oodx2 * (c3inx + dp2dxdy3) \
                    +  hplusy_botleft * oody2 * (c3iny - dp2dxdy3) \
                    +  hplusx_botright * oodx2 * (c4inx + dp2dxdy4) \
                    +  hplusy_botright * oody2 * (c4iny + dp2dxdy4) \
                    +  hcenter_idx * p[idx]

            # if (cnt_x == 1) and (cnt_y == 1):
                # print(lap[idx])
    
            cnt_x += 1
        return lap
    return lap2D_try

def stencil_9pt_3rd_try(elem,node,mpv,ud,diag_inv):
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

    return lambda p : lap2D_3try(p, igx,igy, iicxn, iicyn, hplusx, hplusy, hcenter, oodx2, oody2, x_periodic, y_periodic, x_wall, y_wall, diag_inv)

@jit(nopython=True, nogil=True, cache=True)
def lap2D_3try(p, igx,igy, iicxn, iicyn, hplusx, hplusy, hcenter, oodx2, oody2, x_periodic, y_periodic, x_wall, y_wall, diag_inv):
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


def stencil_9pt_4th_try(elem,node,mpv,ud,diag_inv):
    igx = elem.igx
    igy = elem.igy

    icxn = node.icx
    icyn = node.icy

    iicxn = icxn - (2 * igx)
    iicyn = icyn - (2 * igy)

    iicxn, iicyn = iicyn, iicxn
    ngnc = (iicxn) * (iicyn)

    dx = node.dy
    dy = node.dx

    inner_domain = (slice(igx,-igx),slice(igy,-igy))

    hplusx = mpv.wplus[1][inner_domain]
    hplusy = mpv.wplus[0][inner_domain]

    hcenter = mpv.wcenter[inner_domain]

    diag_inv = diag_inv[inner_domain]

    oodx2 = 0.5 / (dx**2)
    oody2 = 0.5 / (dy**2)
    nine_pt = 0.25 * (2.0) * 1.0

    x_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[0] == BdryType.PERIODIC

    x_wall = ud.bdry_type[1] == BdryType.WALL
    y_wall = ud.bdry_type[0] == BdryType.WALL

    lap = np.zeros((ngnc))

    def lap2D_4try(p):
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

            right = (midright - midmid)
            left = (midmid - midleft)
            top = (midmid - topmid)
            bot = (botmid - midmid)

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
                assert(0)
                
            if y_wall and (cnt_y == (iicyn - 1)):
                hplusx_botleft = 0.
                hplusy_botleft = 0.  
                hplusx_botright = 0.
                hplusy_botright = 0.
                assert(0)
                                 

            dp2dxdy1 = (midmid - midleft) - (topmid - topleft)
            dp2dxdy1 *= nine_pt
            dp2dxdy2 = (midright - midmid) - (topright - topmid)
            dp2dxdy2 *= nine_pt
            dp2dxdy3 = (botmid - botleft) - (midmid - midleft)
            dp2dxdy3 *= nine_pt
            dp2dxdy4 = (botright - botmid) - (midright - midmid)
            dp2dxdy4 *= nine_pt


            lap[idx] = - hplusx_topleft * oodx2 * (left - dp2dxdy1) \
                    -  hplusy_topleft * oody2 * (top - dp2dxdy1) \
                    +  hplusx_topright * oodx2 * (right - dp2dxdy2) \
                    -  hplusy_topright * oody2 * (top + dp2dxdy2) \
                    -  hplusx_botleft * oodx2 * (left + dp2dxdy3) \
                    +  hplusy_botleft * oody2 * (bot - dp2dxdy3) \
                    +  hplusx_botright * oodx2 * (right + dp2dxdy4) \
                    +  hplusy_botright * oody2 * (bot + dp2dxdy4) \
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


            
        return lap
    return lap2D_4try



# class stencil_9pt_5th_try(object):
#     def __init__(self,elem,node,mpv,ud,diag_inv):
#         self.igx = elem.igx
#         self.igy = elem.igy

#         self.icxn = node.icx
#         self.icyn = node.icy

#         self.iicxn = self.icxn - (2 * self.igx)
#         self.iicyn = self.icyn - (2 * self.igy) 

#         self.iicxn , self.iicyn = self.iicyn, self.iicxn 

#         self.ngnc = (self.iicxn) * (self.iicyn)     

#         self.dx = node.dy
#         self.dy = node.dx 

#         inner_domain = (slice(self.igx,-self.igx),slice(self.igy,-self.igy))

#         self.hplusx = mpv.wplus[1][inner_domain].reshape(-1,)
#         self.hplusy = mpv.wplus[0][inner_domain].reshape(-1,)

#         self.hcenter = mpv.wcenter[inner_domain].reshape(-1,)

#         self.diag_inv = diag_inv[inner_domain].reshape(-1,)

#         self.oodx2 = 0.5 / (self.dx**2)
#         self.oody2 = 0.5 / (self.dy**2)
#         self.nine_pt = 0.25 * (2.0) * 1.0

#         self.x_periodic = ud.bdry_type[1] == BdryType.PERIODIC
#         self.y_periodic = ud.bdry_type[0] == BdryType.PERIODIC

#         self.x_wall = ud.bdry_type[1] == BdryType.WALL
#         self.y_wall = ud.bdry_type[0] == BdryType.WALL



#     @staticmethod
#     @jit(nopython=True, nogil=True)
#     def lap2D_3try(p):
#         lap = np.zeros((ngnc))
#         cnt_x = 0
#         cnt_y = 0

#         for idx in range(iicxn * iicyn):
#             ne_topleft = idx - iicxn - 1
#             ne_topright = idx - iicxn 
#             ne_botleft = idx - 1
#             ne_botright = idx

#             # ne_topleft = idx - iicxn
#             # ne_topright = idx - iicxn + 1
#             # ne_botleft = idx
#             # ne_botright = idx +1

#             # ne_topright = idx
#             # ne_topleft = idx - 1
#             # ne_botleft = idx + iicxn - 1
#             # ne_botright = idx + iicxn

#             # ne_topleft = idx
#             # ne_topright = idx + 1
#             # ne_botleft = idx + iicxn
#             # ne_botright = idx + iicxn + 1

#             # get indices of the 9pt stencil
#             topleft_idx = idx - iicxn - 1
#             midleft_idx = idx - 1
#             botleft_idx = idx + iicxn - 1

#             topmid_idx = idx - iicxn
#             midmid_idx = idx
#             botmid_idx = idx + iicxn

#             topright_idx = idx - iicxn + 1
#             midright_idx = idx + 1
#             botright_idx = idx + iicxn + 1

#             if cnt_x == 0:
#                 topleft_idx += iicxn - 1
#                 midleft_idx += iicxn - 1
#                 botleft_idx += iicxn - 1

#                 if x_periodic:
#                     topmid_idx += iicxn - 1
#                     midmid_idx += iicxn - 1
#                     botmid_idx += iicxn - 1

#                 ne_topleft += iicxn - 1
#                 ne_botleft += iicxn - 1

#             if cnt_x == (iicxn - 1):
#                 topright_idx -= iicxn - 1
#                 midright_idx -= iicxn - 1
#                 botright_idx -= iicxn - 1

#                 if x_periodic:
#                     topmid_idx -= iicxn - 1
#                     midmid_idx -= iicxn - 1
#                     botmid_idx -= iicxn - 1

#                 ne_topright -= iicxn - 1
#                 ne_botright -= iicxn - 1

#             if cnt_y == 0:
#                 topleft_idx += ((iicxn) * (iicyn - 1)) 
#                 topmid_idx += ((iicxn) * (iicyn - 1))
#                 topright_idx += ((iicxn) * (iicyn - 1))

#                 if y_periodic:
#                     midleft_idx += ((iicxn) * (iicyn - 1))
#                     midmid_idx += ((iicxn) * (iicyn - 1))
#                     midright_idx += ((iicxn) * (iicyn - 1))

#                 ne_topleft += ((iicxn) * (iicyn - 1))
#                 ne_topright += ((iicxn) * (iicyn - 1))

#             if cnt_y == (iicyn - 1):
#                 botleft_idx -= ((iicxn) * (iicyn - 1))
#                 botmid_idx -= ((iicxn) * (iicyn - 1))
#                 botright_idx -= ((iicxn) * (iicyn - 1))

#                 if y_periodic:
#                     midleft_idx -= ((iicxn) * (iicyn - 1))
#                     midmid_idx -= ((iicxn) * (iicyn - 1))
#                     midright_idx -= ((iicxn) * (iicyn - 1))

#                 ne_botleft -= ((iicxn) * (iicyn - 1))
#                 ne_botright -= ((iicxn) * (iicyn - 1))

#             # if cnt_x == iicxn-1 and cnt_y == 1:
#             #     print([[topleft_idx,topmid,topright],[midleft_idx,midmid,midright],[botleft_idx,botmid,botright]])

#             topleft = p[topleft_idx]
#             midleft = p[midleft_idx]
#             botleft = p[botleft_idx]

#             topmid = p[topmid_idx]
#             midmid = p[midmid_idx]
#             botmid = p[botmid_idx]

#             topright = p[topright_idx]
#             midright = p[midright_idx]
#             botright = p[botright_idx]

#             hplusx_topleft = hplusx[ne_topleft]
#             hplusx_botleft = hplusx[ne_botleft]
#             hplusy_topleft = hplusy[ne_topleft]
#             hplusy_botleft = hplusy[ne_botleft]

#             hplusx_topright = hplusx[ne_topright]
#             hplusx_botright = hplusx[ne_botright]
#             hplusy_topright = hplusy[ne_topright]
#             hplusy_botright = hplusy[ne_botright]

#             right = (midright - midmid)
#             left = (midmid - midleft)
#             top = (midmid - topmid)
#             bot = (botmid - midmid)

#             if x_wall and (cnt_x == 0):
#                 hplusx_topleft = 0.
#                 hplusy_topleft = 0.
#                 hplusx_botleft = 0.
#                 hplusy_botleft = 0.
#             #     # top = 0.
#             #     # bot = 0.
#             #     # right = 0.
#             #     # left= 0.
#             #     # print(hplusx_topright)
#             #     # print(hplusx_botright)

#             if x_wall and (cnt_x == (iicxn - 1)):
#                 hplusx_topright = 0.
#                 hplusy_topright = 0.
#                 hplusx_botright = 0.
#                 hplusy_botright = 0.
#                 # top = 0.
#                 # bot = 0.
#                 # left = 0.
#                 # right = 0.

#             if y_wall and (cnt_y == 0):
#                 hplusx_topleft = 0.
#                 hplusy_topleft = 0. 
#                 hplusx_topright = 0.
#                 hplusy_topright = 0.
#                 assert(0)
                
#             if y_wall and (cnt_y == (iicyn - 1)):
#                 hplusx_botleft = 0.
#                 hplusy_botleft = 0.  
#                 hplusx_botright = 0.
#                 hplusy_botright = 0.
#                 assert(0)
                                 

#             dp2dxdy1 = (midmid - midleft) - (topmid - topleft)
#             dp2dxdy1 *= nine_pt
#             dp2dxdy2 = (midright - midmid) - (topright - topmid)
#             dp2dxdy2 *= nine_pt
#             dp2dxdy3 = (botmid - botleft) - (midmid - midleft)
#             dp2dxdy3 *= nine_pt
#             dp2dxdy4 = (botright - botmid) - (midright - midmid)
#             dp2dxdy4 *= nine_pt


#             lap[idx] = - hplusx_topleft * oodx2 * (left - dp2dxdy1) \
#                     -  hplusy_topleft * oody2 * (top - dp2dxdy1) \
#                     +  hplusx_topright * oodx2 * (right - dp2dxdy2) \
#                     -  hplusy_topright * oody2 * (top + dp2dxdy2) \
#                     -  hplusx_botleft * oodx2 * (left + dp2dxdy3) \
#                     +  hplusy_botleft * oody2 * (bot - dp2dxdy3) \
#                     +  hplusx_botright * oodx2 * (right + dp2dxdy4) \
#                     +  hplusy_botright * oody2 * (bot + dp2dxdy4) \
#                     +  hcenter[idx] * p[idx]
#             # lap[idx] = - hplusx_topleft * oodx2 * ((midmid - midleft) - dp2dxdy1) \
#             #         -  hplusy_topleft * oody2 * ((midmid - topmid) - dp2dxdy1) \
#             #         +  hplusx_topright * oodx2 * ((midright - midmid) - dp2dxdy2) \
#             #         -  hplusy_topright * oody2 * ((midmid - topmid) + dp2dxdy2) \
#             #         -  hplusx_botleft * oodx2 * ((midmid - midleft) + dp2dxdy3) \
#             #         +  hplusy_botleft * oody2 * ((botmid - midmid) - dp2dxdy3) \
#             #         +  hplusx_botright * oodx2 * ((midright - midmid) + dp2dxdy4) \
#             #         +  hplusy_botright * oody2 * ((botmid - midmid) + dp2dxdy4) \
#             #         +  hcenter[idx] * p[idx]

#             # if cnt_x == 0 and x_wall:
#             #     lap[idx] *= 2.
#             # if (cnt_x == iicxn - 1) and x_wall:
#             #     lap[idx] *= 2.

#             # if cnt_y == 0 and y_wall:
#             #     lap[idx] *= 2.
#             # if (cnt_y == iicyn - 1) and y_wall:
#             #     lap[idx] *= 2.

#             lap[idx] *= diag_inv[idx]

#             cnt_x += 1
#             if cnt_x % iicxn == 0:
#                 cnt_y += 1
#                 cnt_x = 0
            
#         return lap
#     return lap2D_5try

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


    