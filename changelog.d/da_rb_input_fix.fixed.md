Fixed the production rising-bubble input (`-ic rb`): its `sol_init` still
called the renamed `hydrostatics.state` (now `integrated_state`), so the case
could not run at all. The MWR-2022 OSSE driver now uses `rb` for the bubble
experiments (native 160x80, t_ref = 1000 s, seeded delth machinery) instead of
scaling up `test_blending_warm_bubble`, which is a 31-step blending smoke and
is CFL-unstable at the paper grid and times. The driver's bubble aux carries
`CFLfixed` but deliberately not `imbal` — the latter now triggers the initial
*hydrostatic* conversion from the hydrostatic-blending work, which NaNs this
nonhydrostatic case; the paper's initial blending is the pseudo-incompressible
one driven by `initial_blending` alone.
