import numba as nb


def get_linop(elem, node, mpv, ud, diag_inv, dt):
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
