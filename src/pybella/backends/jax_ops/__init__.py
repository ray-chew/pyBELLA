"""JAX twins of ``pybella.utils.operators``.

This package is the only place in pyBELLA that imports jax. Equivalence
against the numpy implementation is asserted in float64, so x64 mode is
enabled here, before any jax array can be created.

Contract with the numpy twins:
- same public function names and positional signatures, module-for-module;
- no in-place mutation: functions that mutate arguments in the numpy
  version return the updated arrays instead;
- kernels take flat arrays / floats / static flags (no ModelState objects),
  so they are jit-clean.
"""

from .. import require_jax

require_jax()

import jax

jax.config.update("jax_enable_x64", True)
assert jax.config.jax_enable_x64, "pyBELLA's JAX backend requires float64 (x64) mode"

from . import (
    convolution,
    divergence,
    elliptic_solve,
    finite_difference,
    gradient,
    laplacian,
)
