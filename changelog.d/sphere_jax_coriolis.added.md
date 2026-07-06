Spherical geometry, JAX twins (task 1b, part 1): the hybrid JAX backend
now supports the sphere's Coriolis path. `backends/jax_ops/coriolis.py`
gains a general-H⁻¹ kernel twin (`compute_coefficients_general` /
`apply_inverse_general`) mirroring the numba Sherman-Morrison + (C11)
alpha_w kernel for an arbitrary local up-direction `e`, and now consumes
the spatially varying rotation field (`ud.coriolis_field`) and buoyancy
rank-one term by reusing the numpy `role_components` / `_up_role_components`
(cached, bit-identical). The `coriolis_field` fast-fail guard is removed.
Gated by `test_coriolis_general_sphere_{swe,gw}` in
`test_jax_advection_equivalence.py` (JAX vs numpy to 1e-13 on the TC2
thin-shell and DCMIP-31 shell fixtures).
