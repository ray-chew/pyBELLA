from .. import convolution

def prepare_diag(mpv, node):
    """Highly optimized version with minimal function calls."""
    ndim = node.ndim

    coeff = 0.75 if ndim == 2 else 0.0625 if ndim == 3 else None
    if coeff is None:
        raise ValueError(f"Unsupported ndim: {ndim}")

    dx, dy, dz = node.dx, node.dy, node.dz
    inv_dx2, inv_dy2 = 1.0 / (dx**2), 1.0 / (dy**2)

    diag_kernel = convolution.get_averaging_kernel(ndim, width=2)

    diag = mpv.wcenter.copy()

    # Main diagonal terms
    diag -= (
        coeff * inv_dx2 * convolution.apply_convolution_kernel(mpv.wplus[0], diag_kernel)
    )
    diag -= (
        coeff * inv_dy2 * convolution.apply_convolution_kernel(mpv.wplus[1], diag_kernel)
    )

    if ndim == 2:
        # Cross terms
        inv_dxdy = 1.0 / (dx * dy)
        diag -= (
            coeff
            * inv_dxdy
            * convolution.apply_convolution_kernel(mpv.wplus[0], diag_kernel)
        )
        diag -= (
            coeff
            * inv_dxdy
            * convolution.apply_convolution_kernel(mpv.wplus[1], diag_kernel)
        )
    elif ndim == 3:
        inv_dz2 = 1.0 / (dz**2)
        diag -= (
            coeff
            * inv_dz2
            * convolution.apply_convolution_kernel(mpv.wplus[2], diag_kernel)
        )

    return 1.0 / diag