Terrain phase 4: hydrostates at physical height. `States` gains a field mode
(full per-column fields when terrain is active); `analytical_state` evaluates
its closed form at z(ξ,η) with local dz = J·dη, `integrated_state` gains a
fine-grid quadrature branch. Resting-atmosphere oracle: balanced atmosphere
over a 400 m Agnesi hill stays at rest to ~1e-10 m/s (the solve floor),
identical to flat. `column`/`initial_pressure` assert no terrain.
