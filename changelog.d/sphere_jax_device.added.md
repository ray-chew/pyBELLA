Spherical geometry, JAX twins (task 1, device-resident): the device-
resident backend (`ud.backend="jax-device"`) now runs the sphere entirely
on-device — the fully-fused timestep that actually delivers the GPU/vmap
speedup (the hybrid path round-trips host↔device per kernel). Ported into
`device_kernels.py` / `device_config.py`: the general e_up buoyancy in the
forward step and the alpha_w discard + buoyancy kick in the explicit part;
the general (C11) H^-1 and `ud.coriolis_field` in `_apply_hinv` /
`_coriolis_h_fields` / `_forward_step` (role-ordered field components shipped
in the config); the e_up-parallel stratification coupling in the correction;
the general free-slip WALL mirror + well-balanced gravity fill dispatched in
the device `_ghost_fill` (reusing the shared `BoundaryConfig`); and the
tangent-plane surface constraint at the seven momentum-modifying substeps
(a no-op unless `ud.constrain_to_surface`, so every non-SWE device step
stays bit-identical). The non-vertical-line / `coriolis_field` fast-fail
guards are removed. Gated by `test_window_sphere_{gw,swe}` in
`test_jax_device_equivalence.py` (device vs hybrid to ~1e-16) and
`test_sphere_tc2_device_reproduces_numpy` in `test_jax_device_fullrun.py`
(TC2 device-vs-numpy at the Krylov floor).
