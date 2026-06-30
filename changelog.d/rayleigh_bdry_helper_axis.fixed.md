Fixed a latent axial-agnosticity violation in the Lamb-wave regression cases:
`test_lamb_wave` / `test_unstable_lamb` hard-coded `ud.bdry_type[1]` when applying
the Rayleigh sponge boundary instead of `ud.bdry_type[axes.vertical_axis(ud)]`
(inert today as both run vertical=1, wrong for any other gravity direction). The
four Rayleigh-switch sites now route through a single
`tests/case_setup.apply_rayleigh_bdry` helper (axis from `axes.vertical_axis`),
with `with_tau` folding in the terrain cases' `get_tau_y`. The iny-resizing
`rayleigh_bc_function` stays vertical=1 (documented limitation). Bit-identical at
tol 0 on all four cases.
