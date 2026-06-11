"""Backend selection for the numerical core.

The canonical backend is numpy+numba (``pybella.utils.operators``). The JAX
backend (``pybella.backends.jax_ops``) mirrors it module-for-module and is
validated against it as a golden master; jax is an optional dependency
(``pip install pybella[jax]``) and is only imported inside ``jax_ops``.
"""

import importlib
import importlib.util


def has_jax():
    """True if jax is importable, without importing it."""
    return importlib.util.find_spec("jax") is not None


def require_jax():
    if not has_jax():
        raise ImportError(
            "The JAX backend requires jax; install it with `pip install pybella[jax]`."
        )


def get_operators(backend="numpy"):
    """Return the operators namespace for a backend ("numpy" or "jax").

    The returned module exposes the operator submodules (divergence, gradient,
    convolution, finite_difference, laplacian) with identical public signatures.
    """
    if backend == "numpy":
        # utils.operators has an empty __init__; importing a submodule binds
        # it as an attribute on the parent package.
        for sub in (
            "convolution",
            "divergence",
            "finite_difference",
            "gradient",
            "laplacian.preconditioner",
            "laplacian.lap2D_manual",
            "laplacian.lap3D",
        ):
            importlib.import_module(f"pybella.utils.operators.{sub}")
        return importlib.import_module("pybella.utils.operators")
    if backend == "jax":
        require_jax()
        return importlib.import_module("pybella.backends.jax_ops")
    raise ValueError(f"unknown backend {backend!r}; expected 'numpy' or 'jax'")
