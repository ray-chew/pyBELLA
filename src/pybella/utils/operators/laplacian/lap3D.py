import numpy as np
import numba as nb
from ... import options as opts


def get_linop(elem, node, npf, ud, diag_inv, dt, cij):
    """Build the full-tensor 27-point operator matvec on the node.isc box.

    Discretises  diag_inv * [ sum_ij (1/(d_i d_j)) (1/16) D_i(C_ij F_j p)
    + hcenter p ]  where F_j is the cell-averaged j-derivative, D_i the
    cell-to-node i-difference, and C_ij the (axis-indexed) coefficient
    fields (Gamma^-1 P Theta) * H^-1 — the same H^-1 the momentum
    correction applies, making the elliptic operator consistent with it.
    The legacy operator used bare diagonal coefficients plus hand-coded
    x-z `corrf` cross terms; with H^-1 = identity (no rotation, no
    buoyancy) this operator reproduces it exactly.

    The solve vector is the C-order ravel of an array shaped node.isc
    (interior nodes plus one ghost layer per side, [x, y, z]). The outer
    ghost ring carries zero operator rows; ghost values are reconstructed
    from periodicity inside the kernel.

    cij: 3x3 nested sequence of full cell-shaped coefficient fields.
    """
    oodxyz = node.dxyz
    oodxyz = 1.0 / (oodxyz**2)
    oodx2, oody2, oodz2 = oodxyz[0], oodxyz[1], oodxyz[2]
    odx, ody, odz = 1.0 / node.dx, 1.0 / node.dy, 1.0 / node.dz

    i1 = (slice(1, -1), slice(1, -1), slice(1, -1))

    ndim = elem.ndim
    periodicity = np.empty(ndim, dtype="int64")
    for dim in range(ndim):
        periodicity[dim] = ud.bdry_type[dim] == opts.BdryType.PERIODIC

    # cell-valued coefficient boxes on the (isc - 1) cell box surrounding
    # the node box; copied because the kernel zeroes wall slabs in place
    C = [[np.ascontiguousarray(cij[i][j][i1]) for j in range(3)] for i in range(3)]

    # cross blocks only enter when H^-1 has off-diagonal content
    use_cross = bool(
        max(np.max(np.abs(C[i][j])) for i in range(3) for j in range(3) if i != j) > 0.0
    )

    # unknowns: interior nodes of the box
    hcenter = np.ascontiguousarray(npf.wcenter[i1])
    diag_inv = np.ascontiguousarray(diag_inv)

    # scipy's LinearOperator dtype probe passes an int8 vector; the cast is
    # a no-copy view for the float64 vectors BiCGSTAB actually sends
    return lambda p: lap3D(
        np.asarray(p, dtype=np.float64),
        C[0][0],
        C[0][1],
        C[0][2],
        C[1][0],
        C[1][1],
        C[1][2],
        C[2][0],
        C[2][1],
        C[2][2],
        hcenter,
        oodx2,
        oody2,
        oodz2,
        odx,
        ody,
        odz,
        periodicity,
        diag_inv,
        use_cross,
    )


@nb.jit(nopython=True, cache=False, nogil=False)
def lap3D(
    p0,
    c00,
    c01,
    c02,
    c10,
    c11,
    c12,
    c20,
    c21,
    c22,
    hcenter,
    oodx2,
    oody2,
    oodz2,
    odx,
    ody,
    odz,
    periodicity,
    diag_inv,
    use_cross,
):
    shx, shy, shz = hcenter.shape
    # p0 is the C-order ravel of an [x, y, z] box: reshape must keep that
    # axis order. Copy so the periodic padding below never mutates the
    # caller's (scipy's) vector.
    p = p0.reshape((shx + 2, shy + 2, shz + 2)).copy()

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
            for c in (c00, c01, c02, c10, c11, c12, c20, c21, c22):
                c[0, :, :] = 0.0
                c[-1, :, :] = 0.0
        if bc == True and cnt == 1:
            tmp = p[:, 1, :]
            p[:, 0, :] = p[:, -3, :]
            p[:, -1, :] = p[:, 2, :]
            p[:, 1, :] = p[:, -2, :]
            p[:, -2, :] = tmp
        elif bc == False and cnt == 1:
            for c in (c00, c01, c02, c10, c11, c12, c20, c21, c22):
                c[:, 0, :] = 0.0
                c[:, -1, :] = 0.0
        if bc == True and cnt == 2:
            tmp = p[:, :, 1]
            p[:, :, 0] = p[:, :, -3]
            p[:, :, -1] = p[:, :, 2]
            p[:, :, 1] = p[:, :, -2]
            p[:, :, -2] = tmp
        elif bc == False and cnt == 2:
            for c in (c00, c01, c02, c10, c11, c12, c20, c21, c22):
                c[:, :, 0] = 0.0
                c[:, :, -1] = 0.0
        cnt += 1

    # cell-averaged directional differences F_j(p) on the cell box
    x_fluxes = p[1:, :, :] - p[:-1, :, :]
    y_fluxes = p[:, 1:, :] - p[:, :-1, :]
    z_fluxes = p[:, :, 1:] - p[:, :, :-1]

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

    # diagonal blocks: D_i(C_ii F_i), exactly the legacy structure
    q = c00 * x_flx
    qm = q[:-1, :, :]
    x_flxm = qm[toplefts[0]] + qm[toprights[0]] + qm[botlefts[0]] + qm[botrights[0]]
    qp = q[1:, :, :]
    x_flxp = qp[toplefts[0]] + qp[toprights[0]] + qp[botlefts[0]] + qp[botrights[0]]

    q = c11 * y_flx
    qm = q[:, :-1, :]
    y_flxm = qm[toplefts[1]] + qm[toprights[1]] + qm[botlefts[1]] + qm[botrights[1]]
    qp = q[:, 1:, :]
    y_flxp = qp[toplefts[1]] + qp[toprights[1]] + qp[botlefts[1]] + qp[botrights[1]]

    q = c22 * z_flx
    qm = q[:, :, :-1]
    z_flxm = qm[toplefts[2]] + qm[toprights[2]] + qm[botlefts[2]] + qm[botrights[2]]
    qp = q[:, :, 1:]
    z_flxp = qp[toplefts[2]] + qp[toprights[2]] + qp[botlefts[2]] + qp[botrights[2]]

    lap[1:-1, 1:-1, 1:-1] = (
        oodx2 * coeff * (-x_flxm + x_flxp)
        + oody2 * coeff * (-y_flxm + y_flxp)
        + oodz2 * coeff * (-z_flxm + z_flxp)
        + hcenter * p[1:-1, 1:-1, 1:-1]
    )

    if use_cross:
        cross = np.zeros_like(x_flxm)

        # (i=0, j=1): D_x(C_01 F_y)
        q = c01 * y_flx
        qm = q[:-1, :, :]
        qms = qm[toplefts[0]] + qm[toprights[0]] + qm[botlefts[0]] + qm[botrights[0]]
        qp = q[1:, :, :]
        qps = qp[toplefts[0]] + qp[toprights[0]] + qp[botlefts[0]] + qp[botrights[0]]
        cross += odx * ody * coeff * (qps - qms)

        # (i=0, j=2): D_x(C_02 F_z)
        q = c02 * z_flx
        qm = q[:-1, :, :]
        qms = qm[toplefts[0]] + qm[toprights[0]] + qm[botlefts[0]] + qm[botrights[0]]
        qp = q[1:, :, :]
        qps = qp[toplefts[0]] + qp[toprights[0]] + qp[botlefts[0]] + qp[botrights[0]]
        cross += odx * odz * coeff * (qps - qms)

        # (i=1, j=0): D_y(C_10 F_x)
        q = c10 * x_flx
        qm = q[:, :-1, :]
        qms = qm[toplefts[1]] + qm[toprights[1]] + qm[botlefts[1]] + qm[botrights[1]]
        qp = q[:, 1:, :]
        qps = qp[toplefts[1]] + qp[toprights[1]] + qp[botlefts[1]] + qp[botrights[1]]
        cross += ody * odx * coeff * (qps - qms)

        # (i=1, j=2): D_y(C_12 F_z)
        q = c12 * z_flx
        qm = q[:, :-1, :]
        qms = qm[toplefts[1]] + qm[toprights[1]] + qm[botlefts[1]] + qm[botrights[1]]
        qp = q[:, 1:, :]
        qps = qp[toplefts[1]] + qp[toprights[1]] + qp[botlefts[1]] + qp[botrights[1]]
        cross += ody * odz * coeff * (qps - qms)

        # (i=2, j=0): D_z(C_20 F_x)
        q = c20 * x_flx
        qm = q[:, :, :-1]
        qms = qm[toplefts[2]] + qm[toprights[2]] + qm[botlefts[2]] + qm[botrights[2]]
        qp = q[:, :, 1:]
        qps = qp[toplefts[2]] + qp[toprights[2]] + qp[botlefts[2]] + qp[botrights[2]]
        cross += odz * odx * coeff * (qps - qms)

        # (i=2, j=1): D_z(C_21 F_y)
        q = c21 * y_flx
        qm = q[:, :, :-1]
        qms = qm[toplefts[2]] + qm[toprights[2]] + qm[botlefts[2]] + qm[botrights[2]]
        qp = q[:, :, 1:]
        qps = qp[toplefts[2]] + qp[toprights[2]] + qp[botlefts[2]] + qp[botrights[2]]
        cross += odz * ody * coeff * (qps - qms)

        lap[1:-1, 1:-1, 1:-1] += cross

    lap = lap * diag_inv

    return lap
