NEDAS Phase D4: blending + CFLfixed on the jax-device backend. Blend windows
segment at the conversion steps: the device drivers run the same
`schemes.prepare_blending` as the numpy loop (host round-trip only on
conversion steps), each regime compiles its own static step variant
(`is_compressible` joins the step-cache key), and the psinc variant exports
the predictor half-time pressure the psinc→comp trial extraction reads. The
pinned `initial_blending: True` configs, EnDAB, and the bubble `CFLfixed`
case now run on GPU. Fixed (both sides): the blend trial integration is
capped at exactly one step — previously a last-ULP dt tie could append a
spurious dt≈0 step whose degenerate half-time pressure became the whole
blended dp2n (backend/BLAS-dependent per-member results).
