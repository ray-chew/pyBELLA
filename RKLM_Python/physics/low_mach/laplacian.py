import numpy as np
from inputs.enum_bdry import BdryType
from scipy import sparse
from itertools import product

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

    print(xflux_mid_right.shape)



def stencil_9pt(elem,node,mpv,ud):
    igx = elem.igx
    igy = elem.igy

    icx = elem.icx
    icxn = node.icx
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

    def lap2D(p):
        lap = np.zeros_like(p)

        for j in range(igy - y_periodic, icyn -igy + y_periodic):
            me = j * icx
            mn = j * icxn

            for i in range(igx - x_periodic, icxn -igx + x_periodic):
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
    
