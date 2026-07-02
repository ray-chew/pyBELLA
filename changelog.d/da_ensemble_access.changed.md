Ported the DA layer's member access onto the ModelState API: new
`data_assimilation/ensemble_access.py` is the single place that knows on which
container a DA attribute lives (`CellSolField` vs `NodePressureField`),
replacing the pre-refactor `results[:, loc, ...]` container-index convention
throughout `analysis.py`, `letkf.py`, `etpf.py` and `utils.py`. The `dap.loc`
index dict and magic ghost-pad widths are gone (`elem.igx/igy` now). The
LETKF (Hunt et al. 2007) and ETPF algorithmic cores are verified unchanged
against the pre-refactor reference (`7f0b676~1`) modulo black formatting.
