import numpy as np
from inputs.enum_bdry import BdryType
from scipy import sparse, signal
from numba import jit
import numba as nb

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


def stencil_27pt(elem,node,mpv,ud,diag_inv):
    oodxyz = node.dxyz
    oodxyz = 1./(oodxyz**2)
    oodx2, oody2, oodz2 = oodxyz[0], oodxyz[1], oodxyz[2]

    print(oodxyz)
    i0 = (slice(0,-1),slice(0,-1),slice(0,-1))
    i1 = (slice(1,-1),slice(1,-1),slice(1,-1))
    i2 = (slice(2,-2),slice(2,-2),slice(2,-2))

    ndim = elem.ndim
    periodicity = np.empty(ndim, dtype='int')
    for dim in range(ndim):
        periodicity[dim] = ud.bdry_type[dim] == BdryType.PERIODIC

    # kernel_x = np.ones((1,2,2))
    # kernel_y = np.ones((2,1,2))
    # kernel_z = np.ones((2,2,1))

    hplusx = mpv.wplus[0][i0][i1]
    hplusy = mpv.wplus[1][i0][i1]
    hplusz = mpv.wplus[2][i0][i1]

    # print(hplusx[:,0,:],hplusx[:,-1,:])

    # hplusx = signal.fftconvolve(hplusx,kernel_x,mode='valid')
    # hplusy = signal.fftconvolve(hplusy,kernel_y,mode='valid')
    # hplusz = signal.fftconvolve(hplusz,kernel_z,mode='valid')

    # generalise to support wall along any axis, and not just the y-axis.
    # hplusx[:,0,:],hplusx[:,-1,:] = 0.0, 0.0
    # hplusy[:,0,:],hplusy[:,-1,:] = 0.0, 0.0
    # hplusz[:,0,:],hplusz[:,-1,:] = 0.0, 0.0

    # print(hplusx.shape, hplusy.shape, hplusz.shape)

    hcenter = mpv.wcenter[i2]
    diag_inv = diag_inv[i1]

    return lambda p : lap3D(p, hplusx, hplusy, hplusz, hcenter, oodx2, oody2, oodz2, periodicity, diag_inv)
    # return lambda p : lap3D(p)

@nb.jit(nopython=True, cache=False)
# def lap3D(p0):
def lap3D(p0, hplusx, hplusy, hplusz, hcenter, oodx2, oody2, oodz2, periodicity, diag_inv):
    # p = p0.reshape(shp[0],shp[1],shp[2])
    # p = p0.reshape(hcenter.shape)
    shx, shy, shz = hcenter.shape
    # p = p0.reshape(shx+2,shy+2,shz+2)
    p = p0.reshape(shz+2,shy+2,shx+2)
    # p = p0.reshape(66,34,66)

    coeff = 1./16
    # hplusx, hplusy, hplusz = 1.0, 1.0, 1.0
    # hcenter = 1.0
    # oodx2, oody2, oodz2 = 1.0, 1.0, 1.0
    lap = np.zeros_like(p)

    lefts = [(slice(0,-1),slice(0,None),slice(0,None)),
             (slice(0,None,),slice(0,-1),slice(0,None)),
             (slice(0,None,),slice(0,None),slice(0,-1))
            ]
    
    rights = [(slice(1,None),slice(0,None),slice(0,None)),
              (slice(0,None),slice(1,None),slice(0,None)),
              (slice(0,None),slice(0,None),slice(1,None))
             ]
    
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

    corners = [toplefts,toprights,botlefts,botrights]
    cnt = 0
    for bc in periodicity:
        if bc == True and cnt == 0:
            tmp = p[1,:,:]
            p[0,:,:] = p[-3,:,:]
            p[-1,:,:] = p[2,:,:]
            p[1,:,:] = p[-2,:,:]
            p[-2,:,:] = tmp
        elif bc == True and cnt == 2:
            tmp = p[:,:,1]
            p[:,:,0] = p[:,:,-3]
            p[:,:,-1] = p[:,:,2]
            p[:,:,1] = p[:,:,-2]
            p[:,:,-2] = tmp
        cnt += 1
    
    # pinner = p[1:-1,1:-1,:]
    leftz = p[:,:,:-1]
    rightz = p[:,:,1:]
    
    z_fluxes = (rightz - leftz)

    # pinner = p[1:-1,:,1:-1]
    lefty = p[:,:-1,:]
    righty = p[:,1:,:]

    y_fluxes = (righty - lefty)

    # pinner = p[:,1:-1,1:-1]
    leftx = p[:-1,:,:]
    rightx = p[1:,:,:]

    x_fluxes = (rightx - leftx)

    x_flx = x_fluxes[toplefts[0]] + x_fluxes[toprights[0]] + x_fluxes[botlefts[0]] + x_fluxes[botrights[0]]
    y_flx = y_fluxes[toplefts[1]] + y_fluxes[toprights[1]] + y_fluxes[botlefts[1]] + y_fluxes[botrights[1]]
    z_flx = z_fluxes[toplefts[2]] + z_fluxes[toprights[2]] + z_fluxes[botlefts[2]] + z_fluxes[botrights[2]]

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
                          hcenter * p[1:-1,1:-1,1:-1]

    # lap[1:-1,1:-1,1:-1] = oodx2 * coeff * (-(hplusx[:,:,:-1] * x_fluxes[:,:,:-1]) + (hplusx[:,:,1:] * x_fluxes[:,:,1:])) \
    #     + oody2 * coeff * (-(hplusy[:,:-1,:] * y_fluxes[:,:-1,:]) + (hplusy[:,1:,:] * y_fluxes[:,1:,:])) \
    #     + oodz2 * coeff * (-(hplusx[:-1,:,:] * z_fluxes[:-1,:,:]) + (hplusx[1:,:,:] * z_fluxes[1:,:,:])) \
    #     + hcenter * p[1:-1,1:-1,1:-1]

    lap = lap * diag_inv

    return lap

def precon_diag_prepare(mpv, elem, node, ud):
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

    hplusx = mpv.wplus[0]
    hplusy = mpv.wplus[1]
    if ndim == 3: hplusz = mpv.wplus[2]

    nine_pt = 0.25 * (2.0) * 1.0
    if ndim == 2:
        coeff = 1.0 - nine_pt
    elif ndim == 3:
        coeff = 0.0625
    else:
        assert(0, "ndim = 1?")

    wx = coeff / (dx**2)
    wy = coeff / (dy**2)
    wz = coeff / (dz**2)

    diag_kernel = np.array(np.ones([2] * ndim))

    diag = np.zeros((node.sc)).squeeze()
    diag[idx_n] = -wx * signal.fftconvolve(hplusx[idx_e],diag_kernel,mode='full')[idx_periodic] 
    diag[idx_n] -= wy * signal.fftconvolve(hplusy[idx_e],diag_kernel,mode='full')[idx_periodic]
    if ndim == 3:
        diag[idx_n] -= wz * signal.fftconvolve(hplusz[idx_e],diag_kernel,mode='full')[idx_periodic]

    diag[idx_n] += mpv.wcenter[idx_n]
    diag[idx_n] = 1.0 / diag[idx_n]

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