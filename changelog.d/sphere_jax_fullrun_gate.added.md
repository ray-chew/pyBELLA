Spherical geometry, JAX twins: an in-process Williamson TC2 full-run gate
(`test_sphere_tc2_hybrid_reproduces_numpy` in `test_jax_fullrun.py`)
exercises the whole hybrid-JAX sphere path end to end — general
curvilinear elliptic solve, advection through the metric normals, the
`coriolis_field` H^-1, the general free-slip walls and the tangent-plane
surface constraint. It compares jax-vs-numpy over a short horizon at the
documented initial-projection Krylov floor (~2e-5 in the momenta, above the
1e-5 regression tolerance so the stored-target subprocess gate cannot be
used), far below the ~5e-4 a wrong rotation axis / Coriolis factor gives.
The device-resident ("jax-device") step still fast-fails cleanly on
non-vertical-line metrics (guard moved into `build_device_config` now that
the shared boundary-config no longer rejects them).
