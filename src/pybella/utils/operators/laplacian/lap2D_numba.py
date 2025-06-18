import numpy as np
import numba as nb
from ... import options as opts


def get_linop(npf, node, coriolis, diag_inv, ud):
    dx = node.dx
    dy = node.dy

    hplusx = npf.wplus[0]
    hplusy = npf.wplus[1]
    hcenter = npf.wcenter

    coeffs = [hplusx.T, hplusy.T, hcenter.T]

    shp = node.iisc

    dummy_p = np.zeros((node.isc[1], node.isc[0]))

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
