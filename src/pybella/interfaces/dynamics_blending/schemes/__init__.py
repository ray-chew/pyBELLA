"""Dynamics-regime blending, split from the former single ``schemes.py``.

Public names are re-exported here so existing ``schemes.<name>`` call sites
(``flow_solver/discretisation/time_update.py`` and this package's ``prepare.py``)
are unchanged. Submodules:

* ``blending``       -- the ``Blend`` interface (pressure smoothing, rescale).
* ``comp_psinc``     -- compressible <-> pseudo-incompressible conversions.
* ``swe_lake``       -- shallow-water <-> lake conversions.
* ``hydro_nonhydro`` -- nonhydrostatic <-> hydrostatic conversions; the
  ``do_hydro_to_nonhydro_conv`` body is intentionally kept commented as the
  reinstatement reference (see ``dev_notes/hydrostatic_blending.md``).
* ``orchestration``  -- the per-timestep blending calls invoked by the solver.
"""

from .blending import Blend
from .comp_psinc import do_comp_to_psinc_conv, do_psinc_to_comp_conv
from .swe_lake import do_swe_to_lake_conv, do_lake_to_swe_conv
from .hydro_nonhydro import do_nonhydro_to_hydro_conv, do_hydro_to_nonhydro_conv
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
    "do_nonhydro_to_hydro_conv",
    "do_hydro_to_nonhydro_conv",
    "blending_before_timestep",
    "blending_after_timestep",
    "prepare_blending",
    "check_and_apply_initial_hydrostatic_conversion",
]
