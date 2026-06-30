Routed the last six regression cases (agnesi, Schär, IGW B&B, blending warm
bubble, unstable Lamb, 3D-Coriolis travelling vortex) through
`tests/case_setup.make_diag_state` instead of constructing `DiagnosticState`
inline, so all eleven cases now share the one helper that centralises the
`Nx=inx-1 / Ny=iny-1 / steps=[stepmax-1]` offsets. Case-specific tolerances,
`time_increment`, the f-string `test_name` and the rationale comments are
preserved verbatim as forwarded keywords. Proven value-identical: `vars(diag_state)`
byte-identical for all six, plus bit-identical solver output at tol 0.
