Fixed the psinc->comp blending conversion being silently discarded:
`time_update.do` assigned the pre-conversion `sol`/`npf` aliases returned by
`prepare_blending` back onto the model state, undoing the conversion's rebind
(the pre-refactor code threaded the CONVERTED Sol/mpv through). Also pinned
the throwaway look-ahead steps in `do_psinc_to_comp_conv` /
`do_lake_to_swe_conv` to the limit-regime clock (`window_step = 0`, matching
the paper-era `[0, step-1]` call), so under continuous blending the extracted
half-time pressure is the projected one, not a compressible unprojected one.
The warm-bubble case's loose 1e-0 momenta tolerances — which had hidden the
discarded conversion — are tightened to the 1e-5 defaults, and its golden
master was regenerated deliberately. Visible effect: blended-DA ensembles now
recover the full balanced pressure field after each assimilation (residual
acoustic imbalance in p2_nodes is gone).
