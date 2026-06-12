General curvilinear advection and consumers (tfc pt 4): the advective mass
flux is now uniformly rhoY (N_i . m) / rho on every sweep axis (numpy + JAX
twins share one assembly); CFL uses the contravariant speeds
|N_a . m|/(rho J) and signal bounds c |N_a|/J on all axes; field-mode
hydrostates and ghost-cell hydrostatics use the true vertical thickness
z_eta = J/(N_v)_v instead of J (identical for vertical-line maps). New
end-to-end Tier-2 gates: resting atmosphere + mountain-wave smoke on an
x-stretched general map (test_scripts/test_terrain_stretched_smoke.py).
Legacy G1/G2/z stay stored (consistent with N by construction and asserted
in tests) — deriving them per sweep flip would cost hot-path allocations.
