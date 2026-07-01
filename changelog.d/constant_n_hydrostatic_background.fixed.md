Fixed the hydrostatic background of the two constant-N stratified regression
cases (`test_internal_long_wave`, `test_blending_hydrostatic`): their `sol_init`
now uses the stratification-consistent `hydrostatics.integrated_state` instead of
the ISOTHERMAL `hydrostatics.analytical_state`. The isothermal background gave
these constant-N cases an effective N-squared 3.19x too large (gravity-wave
frequency ~1.79x too fast), a defect masked only because their regression targets
had been self-generated with the same wrong background. The dominant IGW frequency
now matches the thesis-era reference to ~1.3%. Targets for both cases regenerated
(pre-approved bug fix). The five other `alpha_w = 1` gates
(`test_blending_warm_bubble`, `test_travelling_vortex`, `test_igw_baldauf_brdar`,
`test_lamb_wave`, `test_unstable_lamb`) do not use the changed function and were
re-verified bit-identical; `analytical_state` is left untouched (correct for the
genuinely isothermal cases). See `dev_notes/hydrostatic_blending.md` (ROOT CAUSE).
