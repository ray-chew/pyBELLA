"""Dynamics-regime blending, split from the former single ``schemes.py``.

Public names are re-exported here so existing ``schemes.<name>`` call sites
(``flow_solver/discretisation/time_update.py`` and this package's ``prepare.py``)
are unchanged. Submodules:

* ``blending``       -- the ``Blend`` interface (pressure smoothing, rescale).
* ``comp_psinc``     -- compressible <-> pseudo-incompressible conversions.
* ``swe_lake``       -- shallow-water <-> lake conversions.
* ``orchestration``  -- the per-timestep blending calls invoked by the solver.

Nonhydrostatic <-> hydrostatic blending needs NO explicit conversion routine:
the eos schedule (``physics/eos.py``) flips ``is_nonhydrostatic`` and the
hydrostatic elliptic operator + explicit vertical-momentum switch carry the
balance (the thesis-era RKLM_Python mechanism). See
``dev_notes/hydrostatic_blending.md``.
"""

from .blending import Blend
from .comp_psinc import do_comp_to_psinc_conv, do_psinc_to_comp_conv
from .swe_lake import do_swe_to_lake_conv, do_lake_to_swe_conv
from .orchestration import (
    blending_before_timestep,
    blending_after_timestep,
    prepare_blending,
    check_and_apply_initial_hydrostatic_conversion,
)

__all__ = [
    "Blend",
    "do_comp_to_psinc_conv",
    "do_psinc_to_comp_conv",
    "do_swe_to_lake_conv",
    "do_lake_to_swe_conv",
    "blending_before_timestep",
    "blending_after_timestep",
    "prepare_blending",
    "check_and_apply_initial_hydrostatic_conversion",
]
