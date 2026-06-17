"""Small shared helpers for regression-case ``UserData`` construction.

The physics in each ``tests/test_*.py`` ``UserData`` is case-specific and stays
in its own module. These pure helpers factor out only the two genuinely-repeated,
error-prone mechanical idioms so they are written once:

* ``build_bdry`` — the per-instance, object-dtype boundary-type triple (must be a
  fresh array per instance, never a shared class attribute).
* ``make_diag_state`` — the ``DiagnosticState`` wiring, which centralises the
  ``Nx = inx - 1`` / ``Ny = iny - 1`` / ``steps = [stepmax - 1]`` offsets that are
  easy to mistype, while forwarding any case-specific keywords (``plot_compare``,
  ``time_increment``, ``tolerances`` ...).

Neither helper mutates global state; each returns a value the case assigns.
"""

import numpy as np

from ..utils.data_structures import DiagnosticState


def build_bdry(x_bc, y_bc, z_bc):
    """Return a fresh ``(3,)`` object-dtype array of boundary types."""
    bdry = np.empty((3), dtype=object)
    bdry[0] = x_bc
    bdry[1] = y_bc
    bdry[2] = z_bc
    return bdry


def make_diag_state(test_name, file_name, inx, iny, stepmax, steps=None, **kwargs):
    """Build a ``DiagnosticState`` with the standard index/step offsets."""
    return DiagnosticState(
        test_name=test_name,
        file_name=file_name,
        Nx=inx - 1,
        Ny=iny - 1,
        steps=steps if steps is not None else [stepmax - 1],
        **kwargs,
    )
