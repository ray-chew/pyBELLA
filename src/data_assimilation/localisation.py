import numpy as np
from dycore.utils.options import BdryType

def rlocal_5pt(elem,node,ud):
    igx = elem.igx
    igy = elem.igy

    icxn = elem.icx
    icyn = elem.icy

    iicxn = icxn - (2 * igx)
    iicyn = icyn - (2 * igy)

    iicxn, iicyn = iicyn, iicxn

    x_periodic = ud.bdry_type[1] == BdryType.PERIODIC
    y_periodic = ud.bdry_type[0] == BdryType.PERIODIC

    x_wall = ud.bdry_type[1] == BdryType.WALL
    y_wall = ud.bdry_type[0] == BdryType.WALL

    return lambda covar : rlocal_5pt_stencil(covar, iicxn, iicyn, x_periodic, y_periodic, x_wall, y_wall)

# @jit(nopython=True, nogil=True, cache=True)
def rlocal_5pt_stencil(covar, iicxn, iicyn, x_periodic, y_periodic, x_wall, y_wall):
    ngnc = (iicxn) * (iicyn)
    R = np.zeros((ngnc))
    cnt_x = 0
    cnt_y = 0

    for idx in range(iicxn * iicyn):
        # get indices of the 5pt stencil
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

        midleft = covar[midleft_idx]
        topmid = covar[topmid_idx]
        midmid = covar[midmid_idx]
        botmid = covar[botmid_idx]
        midright = covar[midright_idx]

        # if x_wall and (cnt_x == 0):
        #     midleft = 0.0

        # if x_wall and (cnt_x == (iicxn - 1)):
        #     midright = 0.0

        # if y_wall and (cnt_y == 0):
        #     topmid = 0.0
            
        # if y_wall and (cnt_y == (iicyn - 1)):
        #     botmid = 0.0
                                
        R[idx] = (0.5 * midleft + 1.0 * midmid + 0.5 * midright) + (0.5 * topmid + 1.0 * midmid + 0.5 * botmid)

        cnt_x += 1
        if cnt_x % iicxn == 0:
            cnt_y += 1
            cnt_x = 0
        
    return R