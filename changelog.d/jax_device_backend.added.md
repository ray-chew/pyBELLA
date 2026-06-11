Device-resident JAX time loop (`ud.backend = "jax-device"`, env
`PYBELLA_BACKEND=jax-device`): one jitted `step(state, dt)` per time step —
advection sweeps, both explicit/implicit pairs with in-jit bicgstab, every
ghost fill, Rayleigh damping/forcing and diffusion — with host syncs reduced
to one dt scalar per step and full-state pulls at output times (per-step
writer pulls only when a case sets output_timesteps). State is a 7-leaf
pytree; all static structure (geometry, profiles, oriented metric, boundary
and laplacian gather plans) is frozen into a per-window DeviceConfig and
closure-captured; exactly two compiled step variants (Strang parity), cached
across output windows. Unsupported configs (blending, ArakawaKonor, acoustic
dt, file forcing, debug writers) raise with an actionable message. All 10
regression cases pass against the stored golden masters; device-vs-hybrid
window equivalence at the documented Krylov floor; ~3.5x faster than numpy
per step already at 64^2 on CPU. Also: `run_scripts/bench_device.py`
benchmark harness and a `pybella[jax-cuda]` install extra for GPU clusters
(code is hardware-agnostic).
