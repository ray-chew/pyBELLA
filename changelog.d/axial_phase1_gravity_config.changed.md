Axial-agnosticity Phase 1: `gravity_direction` is now a configuration input
(default 1 = y-vertical, validated; 2D runs require 1), `gravity_strength` and
`coriolis_strength` are computed role-based from it, and all seven hardcoded
`ud.gravity_strength[1]` reads route through `axes.vertical_axis(ud)`.
Bit-identical for all existing cases (verified at tol=0 on full run outputs).
