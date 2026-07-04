# NEDAS-vs-native comparison figures

Scripts that generated the cross-engine comparison figures (member contour
panels, imshow variants, RMSE/spread overlays incl. the t=0 initial error).
Inputs: the NEDAS memory checkpoints under `outputs/nedas_*/` (worktree) and
the native ensemble/truth H5s + `osse_final_*` diagnostic CSVs (produced by
`run_scripts/osse_mwr2022.py` + `osse_diagnostics.py`). All regenerable from
seeds; see dev_notes/nedas_interface.md for the run recipes.

- plot_nedas_vs_native_contours.py — TV p2 member contours (momentum + all-quantities)
- plot_all_imshow.py               — TV all-quantities, imshow variant
- plot_bubble_members.py           — RB p2 member contours at t=1.0
- plot_rmse_overlays.py            — per-field RMSE/spread overlays, native vs NEDAS, incl. t=0
