Lifted the JAX boundary's field-mode + incompressible gravity-fill guard for
the general (spherical) metric, so the TC2-style initial projection is now
available on the JAX device/hybrid backends for the deep compressible shell
(Hughes & Jablonowski). `jax_ops/boundary._general_gravity_ops` now gathers
the field-mode `HydroState.rhoY0` at the image cell as a jnp field slice — the
twin of numpy's `cell_boundary._hydro_at` — oriented with a plain transpose
(`_orient_leaf`) since the hydrostate is never sweep-flipped. Because the
gravity fill bakes the `is_compressible` regime in statically and
`do_initial_projection` toggles it on the same `ud`, `get_boundary_config` now
also invalidates on a compressibility change (`BoundaryConfig.compressible_ref`)
so a projection-era incompressible config is not reused for the compressible
device loop. Validated bit-exact against numpy on the same input
(`test_field_mode_gravity_ghost_device_matches_numpy`, ghost diff 2e-16); the
sphere TC2 hybrid/device gates, JAX boundary-equivalence suite, and the
travelling-vortex (projection) + agnesi (terrain) device fullruns against
golden targets all still pass. The vertical-line terrain twin (`_gravity_ops`)
keeps the guard — no vertical-line case runs an incompressible fill. NOTE:
end-to-end projection on the deep compressible H&J shell is ill-conditioned
(its bicgstab floor makes the backends diverge ~10% in the momenta); the fill
is exact, and the case runs projection-free by default regardless.
