"""Terrain h == 0 identity oracle.

Two guarantees, two mechanisms:

1. **Bypass identity (bit-exact).** No registered golden-master case
   defines ``orography``, so ``elem.metric is None`` and every solver
   call site takes the pre-terrain code path untouched. The 9-case
   golden-master suite (``test_flow_solver.py``) is the authority; this
   file just pins the precondition.

2. **Forced-flat identity (~1e-13).** The metric machinery switched ON
   with ``h == 0`` is algebraically the identity (J == 1, G == 0) but
   multiplies extra factors through the operators. As each metric-aware
   operator lands (Phases 1-6 of the terrain plan), a comparison of the
   flat-metric path against the plain path is added here — run in both
   the compressible and pseudo-incompressible regimes to cover the
   pi-update / wcenter consistency factors.
"""

import importlib

import numpy as np

from pybella.interfaces.ic_config import IC_MODULES
from pybella.utils import user_data

# the golden-master regression cases (smoke_* cases are exempt: they own
# no target and are allowed to grow terrain)
GOLDEN_MASTER_CASES = [
    "test_travelling_vortex",
    "test_travelling_vortex_3d_coriolis",
    "test_internal_long_wave",
    "test_igw_baldauf_brdar",
    "test_lamb_wave",
    "test_blending_warm_bubble",
    "test_unstable_lamb",
    "test_swe_vortex",
    "test_straka",
]


def test_golden_master_cases_have_no_orography():
    """Precondition for bit-identity: every target-bearing case must take
    the metric-bypass path (elem.metric is None)."""
    for ic in GOLDEN_MASTER_CASES:
        module = importlib.import_module(IC_MODULES[ic])
        ud = user_data.UserDataInit(**vars(module.UserData()))
        assert getattr(ud, "orography", None) is None, (
            f"{ic} defines orography — golden-master cases must stay on "
            "the uniform-Cartesian bypass path"
        )


# --- forced-flat operator comparisons (fleshed out per terrain phase) -------
#
# Phase 1: divergence (J grad F vs plain kernel, flat metric)
# Phase 2: pressure gradients (A-map applied with J=1, G=0)
# Phase 3: elliptic operator C_ij and wcenter
# Phase 6: advective fluxes + CFL
