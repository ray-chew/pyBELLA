"""PyBellaObs — synthetic obs for the MWR-2022 parity OSSEs (N1 skeleton).

Two modes (dataset_def.pybella keys), see dev_notes/nedas_interface.md:

1. obs_file mode (default for N2 parity): read the NATIVE obs HDF5 produced by
   run_scripts/osse_mwr2022.py (labels ``<attr>_ensemble_mem=0_<t:.3f>_
   after_full_step``, times rounded to 3 decimals) and expose the values at
   'prescribed' inner-cell-centre positions, so NEDAS and the native LETKF
   assimilate byte-identical observations.
2. regenerate mode: reproduce the native design in-framework — obs at inner
   cell centres, every-10th-point sparsity (seed 777), Gaussian noise
   (seed 888) with err std = 5% field variance (VarCov).
"""

import numpy as np
from NEDAS.datasets.synthetic import SyntheticObs


class PyBellaObs(SyntheticObs):
    obs_file: str | None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # TODO(N1): VarDesc entries for 'momentum' (vector) mirroring
        # PyBellaModel.variables; obs values are point samples of the state,
        # so no custom obs_operator is needed (identity/interp default).

    def generate_obs_network(self, **kwargs):
        """Return obs_seq dict (obs/t/z/y/x/err_std) per the modes above."""
        raise NotImplementedError
