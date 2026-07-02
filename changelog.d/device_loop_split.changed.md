Split the jax-device host window driver (`run_window`, CFL/dt host control,
support guard, forcing eval, compile cache) out of `device_kernels.py` into
`backends/jax_ops/device_loop.py`; `device_step` re-exports are unchanged.
