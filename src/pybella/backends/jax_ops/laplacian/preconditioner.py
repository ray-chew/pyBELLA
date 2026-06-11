"""JAX twin of :mod:`pybella.utils.operators.laplacian.preconditioner`.

Same signature as the numpy twin: ``npf``/``node`` are read-only here and
this runs once per linear-system setup (host-side construction time), so
taking the objects is jit-safe.
"""

import jax.numpy as jnp

from .. import convolution


def prepare_diag(npf, node, cii=None):
    """Diagonal preconditioner for the Helmholtz solve.

    cii: optional per-axis diagonal coefficient fields (the C_ii of the
    full-tensor 3D operator); defaults to npf.wplus.
    """
    ndim = node.ndim
    w0, w1, w2 = (
        (npf.wplus[0], npf.wplus[1], npf.wplus[2] if ndim == 3 else None)
        if cii is None
        else cii
    )

    coeff = 0.75 if ndim == 2 else 0.0625 if ndim == 3 else None
    if coeff is None:
        raise ValueError(f"Unsupported ndim: {ndim}")

    dx, dy, dz = node.dx, node.dy, node.dz
    inv_dx2, inv_dy2 = 1.0 / (dx**2), 1.0 / (dy**2)

    diag_kernel = convolution.get_averaging_kernel(ndim, width=2)

    diag = jnp.asarray(npf.wcenter)

    diag = diag - coeff * inv_dx2 * convolution.apply_convolution_kernel(
        w0, diag_kernel
    )
    diag = diag - coeff * inv_dy2 * convolution.apply_convolution_kernel(
        w1, diag_kernel
    )

    if ndim == 2:
        inv_dxdy = 1.0 / (dx * dy)
        diag = diag - coeff * inv_dxdy * convolution.apply_convolution_kernel(
            w0, diag_kernel
        )
        diag = diag - coeff * inv_dxdy * convolution.apply_convolution_kernel(
            w1, diag_kernel
        )
    elif ndim == 3:
        inv_dz2 = 1.0 / (dz**2)
        diag = diag - coeff * inv_dz2 * convolution.apply_convolution_kernel(
            w2, diag_kernel
        )

    return 1.0 / diag
