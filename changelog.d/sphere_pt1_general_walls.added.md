Spherical geometry, stage A2: general curvilinear free-slip walls — the
ghost-cell fill mirrors the CONTRAVARIANT momentum triple with the local
area normals (per-cell 3x3 solve; exact wall-flux cancellation) for
non-vertical-line metrics, in both the no-gravity WALL path and the
gravity path (tangential-velocity copy + `e_up` normal reflection,
`h_v`-based ghost spacing). Vertical-line maps keep the legacy slope-term
recipe verbatim (bit-identity); the JAX ghost fill fast-fails on
spherical metrics until the SWE stage.
