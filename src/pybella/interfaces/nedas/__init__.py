"""NEDAS adapter for pyBELLA.

NEDAS (github.com/nansencenter/NEDAS, pinned ==1.2.0) drives pyBELLA as an
ensemble-DA "model": per-member ModelState forecasts via
flow_solver.discretisation.time_update.do, in-memory state exchange
(io_mode: online), and its batch (LETKF-like) / ETKF assimilators.

NEDAS is an optional extra (``pip install "pybella[nedas]"``); nothing in the
deterministic solver imports this package. NEDAS v1.2.0 resolves model and
dataset classes from hardcoded registries that import ``NEDAS.models.<name>``,
so registration injects our modules into ``sys.modules`` under those names —
call :func:`register` before constructing a NEDAS ``Context``/scheme.
"""


def register() -> None:
    """Register PyBellaModel and PyBellaObs with the NEDAS v1.2.0 registries."""
    import sys

    import NEDAS.datasets
    import NEDAS.models

    from . import model as _model_module
    from . import obs as _obs_module

    NEDAS.models.registry["pybella"] = "PyBellaModel"
    sys.modules["NEDAS.models.pybella"] = _model_module

    NEDAS.datasets.registry["pybella"] = "PyBellaObs"
    sys.modules["NEDAS.datasets.pybella"] = _obs_module


__all__ = ["register"]
