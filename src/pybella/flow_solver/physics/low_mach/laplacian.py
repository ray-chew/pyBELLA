import numpy as np
import scipy as sp
import numba as nb

from ....utils import options as opts
from ....utils import operators


def get_lap2D_stencil(mpv, node, coriolis, diag_inv, ud):
    dx = node.dx
    dy = node.dy

    hplusx = mpv.wplus[0]
    hplusy = mpv.wplus[1]
    hcenter = mpv.wcenter

    coeffs = [hplusx.T, hplusy.T, hcenter.T]

    shp = node.iisc

    dummy_p = np.zeros((node.isc[1], node.isc[0]))

    # if hasattr(ud, "ATMOSPHERIC_EXTENSION"):
    #     return lambda p: lap2D_extended(
    #         p, dummy_p, dx, dy, coeffs, diag_inv.T, coriolis, shp
    #     )

    if hasattr(ud, "ATMOSPHERIC_EXTENSION") and ud.ATMOSPHERIC_EXTENSION:
        boundary_handler = periodic_x_wall_y
    else:
        boundary_handler = periodic

        x_wall = (
            ud.bdry_type[0] == opts.BdryType.WALL
            or ud.bdry_type[0] == opts.BdryType.RAYLEIGH
        )
        y_wall = (
            ud.bdry_type[1] == opts.BdryType.WALL
            or ud.bdry_type[1] == opts.BdryType.RAYLEIGH
        )

        if x_wall:
            coeffs[0], coeffs[1] = apply_x_wall_boundary_coeffs(coeffs[0], coeffs[1])
        if y_wall:
            coeffs[0], coeffs[1] = apply_y_wall_boundary_coeffs(coeffs[0], coeffs[1])

    return lambda p: lap2D_generic(
        p, dummy_p, dx, dy, coeffs, diag_inv.T, coriolis, shp, boundary_handler
    )

@nb.njit(cache=True)
def lap2D_generic(p, dp, dx, dy, coeffs, diag_inv, coriolis, shp, boundary_handler):
    p = p.reshape(shp[1], shp[0])
    dp[1:-1, 1:-1] = p

    dp[...] = boundary_handler(dp)
    dp[...] = kernel_9pt(
        dp,
        dx,
        dy,
        coeffs[0],
        coeffs[1],
        coeffs[2],
        diag_inv,
        coriolis[0],
        coriolis[1],
        coriolis[2],
        coriolis[3],
    )

    return dp[1:-1, 1:-1].ravel()


@nb.njit(cache=True)
def periodic_x_wall_y(arr):
    # periodic padding
    arr[:, 0] = arr[:, -3]
    arr[:, -1] = arr[:, 2]

    # wall padding
    arr[0, :] = arr[2, :]
    arr[-1, :] = arr[-3, :]

    return arr

@nb.njit(cache=True)
def periodic(arr):
    """Apply periodic boundary conditions"""
    arr[:, 0] = arr[:, -3]
    arr[:, -1] = arr[:, 2]
    arr[0, :] = arr[-3, :]
    arr[-1, :] = arr[2, :]

    return arr

@nb.njit(cache=True)
def apply_x_wall_boundary_coeffs(hpx, hpy):
    """Apply wall boundary conditions by modifying coefficients"""
    hpx[:, 0] = 0.0
    hpy[:, 0] = 0.0
    hpx[:, -1] = 0.0
    hpy[:, -1] = 0.0
    hpx[:, 1] = 0.0
    hpy[:, 1] = 0.0
    hpx[:, -2] = 0.0
    hpy[:, -2] = 0.0

    return hpx, hpy


@nb.njit(cache=True)
def apply_y_wall_boundary_coeffs(hpx, hpy):
    """Apply wall boundary conditions by modifying coefficients"""
    hpx[0, :] = 0.0
    hpy[0, :] = 0.0
    hpx[-1, :] = 0.0
    hpy[-1, :] = 0.0
    hpx[1, :] = 0.0
    hpy[1, :] = 0.0
    hpx[-2, :] = 0.0
    hpy[-2, :] = 0.0

    return hpx, hpy

@nb.stencil
def kernel_9pt(a, dx, dy, hpx, hpy, hpc, diag_inv, cxx, cyy, cxy, cyx):
    oodx = 1.0 / dx
    oody = 1.0 / dy

    topleft = a[1, -1]
    topmid = a[1, 0]
    topright = a[1, 1]

    midleft = a[0, -1]
    midmid = a[0, 0]
    midright = a[0, 1]

    botleft = a[-1, -1]
    botmid = a[-1, 0]
    botright = a[-1, 1]

    hpx_bl = hpx[0, 0]
    hpx_br = hpx[0, 1]
    hpx_tl = hpx[1, 0]
    hpx_tr = hpx[1, 1]

    hpy_bl = hpy[0, 0]
    hpy_br = hpy[0, 1]
    hpy_tl = hpy[1, 0]
    hpy_tr = hpy[1, 1]

    cxx_bl = cxx[0, 0]
    cxx_br = cxx[0, 1]
    cxx_tl = cxx[1, 0]
    cxx_tr = cxx[1, 1]

    cyy_bl = cyy[0, 0]
    cyy_br = cyy[0, 1]
    cyy_tl = cyy[1, 0]
    cyy_tr = cyy[1, 1]

    cxy_bl = cxy[0, 0]
    cxy_br = cxy[0, 1]
    cxy_tl = cxy[1, 0]
    cxy_tr = cxy[1, 1]

    cyx_bl = cyx[0, 0]
    cyx_br = cyx[0, 1]
    cyx_tl = cyx[1, 0]
    cyx_tr = cyx[1, 1]

    Dx_tl = 0.5 * (topmid - topleft + midmid - midleft) * hpx_tl
    Dx_tr = 0.5 * (topright - topmid + midright - midmid) * hpx_tr
    Dx_bl = 0.5 * (botmid - botleft + midmid - midleft) * hpx_bl
    Dx_br = 0.5 * (botright - botmid + midright - midmid) * hpx_br

    Dy_tl = 0.5 * (topmid - midmid + topleft - midleft) * hpy_tl
    Dy_tr = 0.5 * (topright - midright + topmid - midmid) * hpy_tr
    Dy_bl = 0.5 * (midmid - botmid + midleft - botleft) * hpy_bl
    Dy_br = 0.5 * (midright - botright + midmid - botmid) * hpy_br

    Dxx = (
        0.5
        * (cxx_tr * Dx_tr - cxx_tl * Dx_tl + cxx_br * Dx_br - cxx_bl * Dx_bl)
        * oodx
        * oodx
    )
    Dyy = (
        0.5
        * (cyy_tr * Dy_tr - cyy_br * Dy_br + cyy_tl * Dy_tl - cyy_bl * Dy_bl)
        * oody
        * oody
    )
    Dyx = (
        0.5
        * (cxy_br * Dy_br - cxy_bl * Dy_bl + cxy_tr * Dy_tr - cxy_tl * Dy_tl)
        * oody
        * oodx
    )
    Dxy = (
        0.5
        * (cyx_tr * Dx_tr - cyx_br * Dx_br + cyx_tl * Dx_tl - cyx_bl * Dx_bl)
        * oodx
        * oody
    )

    return ((Dxx + Dyy + Dyx + Dxy) + hpc[0, 0] * a[0, 0]) * diag_inv[0, 0]

def stencil_27pt(elem, node, mpv, ud, diag_inv, dt):
    oodxyz = node.dxyz
    oodxyz = 1.0 / (oodxyz**2)
    oodx2, oody2, oodz2 = oodxyz[0], oodxyz[1], oodxyz[2]
    odx, odz = 1.0 / node.dx, 1.0 / node.dz

    i0 = (slice(0, -1), slice(0, -1), slice(0, -1))
    i1 = (slice(1, -1), slice(1, -1), slice(1, -1))
    i2 = (slice(2, -2), slice(2, -2), slice(2, -2))

    ndim = elem.ndim
    periodicity = np.empty(ndim, dtype="int")
    for dim in range(ndim):
        periodicity[dim] = ud.bdry_type[dim] == opts.BdryType.PERIODIC

    hplusx = mpv.wplus[0][i0][i1]
    hplusy = mpv.wplus[1][i0][i1]
    hplusz = mpv.wplus[2][i0][i1]

    hcenter = mpv.wcenter[i2]
    diag_inv = diag_inv[i1]

    corrf = dt * ud.coriolis_strength[0]

    return lambda p: lap3D(
        p,
        hplusx,
        hplusy,
        hplusz,
        hcenter,
        oodx2,
        oody2,
        oodz2,
        periodicity,
        diag_inv,
        corrf,
        odx,
        odz,
    )


@nb.jit(nopython=True, cache=False, nogil=False)
def lap3D(
    p0,
    hplusx,
    hplusy,
    hplusz,
    hcenter,
    oodx2,
    oody2,
    oodz2,
    periodicity,
    diag_inv,
    corrf,
    odx,
    odz,
):
    shx, shy, shz = hcenter.shape
    p = p0.reshape(shz + 2, shy + 2, shx + 2)

    coeff = 1.0 / 16
    lap = np.zeros_like(p)

    # cut out four cubes from the 3d array corresponding to the nodes... in each axial direction.
    toplefts = [
        (slice(0, None), slice(0, -1), slice(0, -1)),
        (slice(0, -1), slice(0, None), slice(0, -1)),
        (slice(0, -1), slice(0, -1), slice(0, None)),
    ]
    toprights = [
        (slice(0, None), slice(0, -1), slice(1, None)),
        (slice(1, None), slice(0, None), slice(0, -1)),
        (slice(0, -1), slice(1, None), slice(0, None)),
    ]

    botlefts = [
        (slice(0, None), slice(1, None), slice(0, -1)),
        (slice(0, -1), slice(0, None), slice(1, None)),
        (slice(1, None), slice(0, -1), slice(0, None)),
    ]
    botrights = [
        (slice(0, None), slice(1, None), slice(1, None)),
        (slice(1, None), slice(0, None), slice(1, None)),
        (slice(1, None), slice(1, None), slice(0, None)),
    ]

    cnt = 0
    for bc in periodicity:
        if bc == True and cnt == 0:
            tmp = p[1, :, :]
            p[0, :, :] = p[-3, :, :]
            p[-1, :, :] = p[2, :, :]
            p[1, :, :] = p[-2, :, :]
            p[-2, :, :] = tmp
        elif bc == False and cnt == 0:
            hplusx[0, :, :] = 0.0
            hplusx[-1, :, :] = 0.0
            hplusy[0, :, :] = 0.0
            hplusy[-1, :, :] = 0.0
            hplusz[0, :, :] = 0.0
            hplusz[-1, :, :] = 0.0
        if bc == True and cnt == 1:
            tmp = p[:, 1, :]
            p[:, 0, :] = p[:, -3, :]
            p[:, -1, :] = p[:, 2, :]
            p[:, 1, :] = p[:, -2, :]
            p[:, -2, :] = tmp
        elif bc == False and cnt == 1:
            hplusx[:, 0, :] = 0.0
            hplusx[:, -1, :] = 0.0
            hplusy[:, 0, :] = 0.0
            hplusy[:, -1, :] = 0.0
            hplusz[:, 0, :] = 0.0
            hplusz[:, -1, :] = 0.0
        if bc == True and cnt == 2:
            tmp = p[:, :, 1]
            p[:, :, 0] = p[:, :, -3]
            p[:, :, -1] = p[:, :, 2]
            p[:, :, 1] = p[:, :, -2]
            p[:, :, -2] = tmp
        elif bc == False and cnt == 2:
            hplusx[:, :, 0] = 0.0
            hplusx[:, :, -1] = 0.0
            hplusy[:, :, 0] = 0.0
            hplusy[:, :, -1] = 0.0
            hplusz[:, :, 0] = 0.0
            hplusz[:, :, -1] = 0.0
        cnt += 1

    leftz = p[:, :, :-1]
    rightz = p[:, :, 1:]

    z_fluxes = rightz - leftz

    lefty = p[:, :-1, :]
    righty = p[:, 1:, :]

    y_fluxes = righty - lefty

    leftx = p[:-1, :, :]
    rightx = p[1:, :, :]

    x_fluxes = rightx - leftx

    x_flx = (
        x_fluxes[toplefts[0]]
        + x_fluxes[toprights[0]]
        + x_fluxes[botlefts[0]]
        + x_fluxes[botrights[0]]
    )
    y_flx = (
        y_fluxes[toplefts[1]]
        + y_fluxes[toprights[1]]
        + y_fluxes[botlefts[1]]
        + y_fluxes[botrights[1]]
    )
    z_flx = (
        z_fluxes[toplefts[2]]
        + z_fluxes[toprights[2]]
        + z_fluxes[botlefts[2]]
        + z_fluxes[botrights[2]]
    )

    hxzp = hplusx * z_flx
    hxzpm = hxzp[:-1, :, :]
    hxzpm = (
        hxzpm[toplefts[0]]
        + hxzpm[toprights[0]]
        + hxzpm[botlefts[0]]
        + hxzpm[botrights[0]]
    )
    hxzpp = hxzp[1:, :, :]
    hxzpp = (
        hxzpp[toplefts[0]]
        + hxzpp[toprights[0]]
        + hxzpp[botlefts[0]]
        + hxzpp[botrights[0]]
    )

    hzxp = hplusz * x_flx
    hzxpm = hzxp[:, :, :-1]
    hzxpm = (
        hzxpm[toplefts[2]]
        + hzxpm[toprights[2]]
        + hzxpm[botlefts[2]]
        + hzxpm[botrights[2]]
    )
    hzxpp = hzxp[:, :, 1:]
    hzxpp = (
        hzxpp[toplefts[2]]
        + hzxpp[toprights[2]]
        + hzxpp[botlefts[2]]
        + hzxpp[botrights[2]]
    )

    x_flx = hplusx * x_flx
    x_flxm = x_flx[:-1, :, :]
    x_flxm = (
        x_flxm[toplefts[0]]
        + x_flxm[toprights[0]]
        + x_flxm[botlefts[0]]
        + x_flxm[botrights[0]]
    )
    x_flxp = x_flx[1:, :, :]
    x_flxp = (
        x_flxp[toplefts[0]]
        + x_flxp[toprights[0]]
        + x_flxp[botlefts[0]]
        + x_flxp[botrights[0]]
    )

    y_flx = hplusy * y_flx
    y_flxm = y_flx[:, :-1, :]
    y_flxm = (
        y_flxm[toplefts[1]]
        + y_flxm[toprights[1]]
        + y_flxm[botlefts[1]]
        + y_flxm[botrights[1]]
    )
    y_flxp = y_flx[:, 1:, :]
    y_flxp = (
        y_flxp[toplefts[1]]
        + y_flxp[toprights[1]]
        + y_flxp[botlefts[1]]
        + y_flxp[botrights[1]]
    )

    z_flx = hplusz * z_flx
    z_flxm = z_flx[:, :, :-1]
    z_flxm = (
        z_flxm[toplefts[2]]
        + z_flxm[toprights[2]]
        + z_flxm[botlefts[2]]
        + z_flxm[botrights[2]]
    )
    z_flxp = z_flx[:, :, 1:]
    z_flxp = (
        z_flxp[toplefts[2]]
        + z_flxp[toprights[2]]
        + z_flxp[botlefts[2]]
        + z_flxp[botrights[2]]
    )

    lap[1:-1, 1:-1, 1:-1] = (
        oodx2 * coeff * (-x_flxm + x_flxp)
        + oody2 * coeff * (-y_flxm + y_flxp)
        + oodz2 * coeff * (-z_flxm + z_flxp)
        + +1.0 * odx * odz * coeff * corrf * (hxzpp - hxzpm)
        + -1.0 * odx * odz * coeff * corrf * (hzxpp - hzxpm)
        + hcenter * p[1:-1, 1:-1, 1:-1]
    )

    lap = lap * diag_inv

    return lap

def precon_diag_prepare(mpv, node):
    """Highly optimized version with minimal function calls."""
    ndim = node.ndim
    
    coeff = 0.75 if ndim == 2 else 0.0625 if ndim == 3 else None
    if coeff is None:
        raise ValueError(f"Unsupported ndim: {ndim}")
    
    dx, dy, dz = node.dx, node.dy, node.dz
    inv_dx2, inv_dy2 = 1.0 / (dx**2), 1.0 / (dy**2)
    
    diag_kernel = operators.get_averaging_kernel(ndim, width=2)

    diag = mpv.wcenter.copy()
    
    # Main diagonal terms
    diag -= coeff * inv_dx2 * operators.apply_convolution_kernel(mpv.wplus[0], diag_kernel)
    diag -= coeff * inv_dy2 * operators.apply_convolution_kernel(mpv.wplus[1], diag_kernel)
    
    if ndim == 2:
        # Cross terms
        inv_dxdy = 1.0 / (dx * dy)
        diag -= coeff * inv_dxdy * operators.apply_convolution_kernel(mpv.wplus[0], diag_kernel)
        diag -= coeff * inv_dxdy * operators.apply_convolution_kernel(mpv.wplus[1], diag_kernel)
    elif ndim == 3:
        inv_dz2 = 1.0 / (dz**2)
        diag -= coeff * inv_dz2 * operators.apply_convolution_kernel(mpv.wplus[2], diag_kernel)
    
    return 1.0 / diag
