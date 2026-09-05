0.51.0 (2026-09-04)
-------------------

Infrastructure
^^^^^^^^^^^^^^

- Added `test_scripts/test_da_smoke.py` to CI: a seeded 2-cycle, N=4
  travelling-vortex OSSE for both LETKF (rloc) and ETPF, asserting that the
  analysis ensemble mean beats the forecast against the regenerated truth on
  the observed momentum fields and that the analysis spread contracts.
  Statistical, not bitwise (the dask-chunked rloc analysis makes bitwise
  comparisons fragile). The DA layer is frozen at the MWR-2022 reproduction:
  2D x-y, vertical=1, numpy backend, bug fixes only — superseded by the NEDAS
  interface as the maintained engine. (da_ci_smoke)
- Repro gate: new `blending` case set (warm bubble + hydrostatic) and
  `test_blending_hydrostatic` added to the `full` set. (repro_gate_blending_set)
- Added ``test_scripts/repro_gate.py``, a reproducibility-gate harness for
  behaviour-preserving refactors: ``capture`` snapshots each case's output H5 and
  per-field max-abs to a scratch baseline, ``check`` reruns and requires the inline
  CompareSol gate to stay green *and* the output to be bit-identical (tol 0) to the
  baseline. Bit-identity is the strong oracle proving a refactor changed no numerics. (repro_gate_harness)
- Physics-vs-reference golden masters for the pole-to-pole sphere, plus the CI
  wiring that was missing for the Stage F (poles) and TFC-generalization gates.

  `test_sphere_swe_tc2_global` and `test_sphere_swe_tc6_global` now have stored
  references (`outputs/target_sphere_swe_tc2_global/`,
  `outputs/target_sphere_swe_tc6_global/`, 128 x 72, step 30) and are registered
  in `test_scripts/test_flow_solver.py`. Until now the only CI coverage of the
  global sphere path was `test_jax_fullrun` / `test_jax_device_fullrun`, which
  gate jax-vs-numpy EQUIVALENCE — a numpy pole regression moved both backends
  together and passed silently. Fresh-target verification is clean at ~1e-8
  (rho 3.0e-8, momenta 2.1e-8, rhoY 7.2e-8, p2_nodes 7.2e-8 for TC2-global;
  ~1.1e-7 worst for TC6-global) against the 1e-5 gate, i.e. ~100x headroom over
  the run-to-run bicgstab scatter.

  Fixed a latent naming bug that made the comparator unusable on all four
  `*_global` cases: their `DiagnosticState` was built from `(inx, iny)` while
  `output_suffix` uses `(inx, inz)`, so `CompareSol` looked for a reference
  filename the writer never emits (`..._128_1_stripped.h5` vs the written
  `..._128_72_stripped.h5`). The two SWE global cases now pass `inz`, and both
  also carry their own `diag_updt_targets` switch so regenerating a channel
  target no longer regenerates the global one (and vice versa).

  CI: the integration job additionally runs the sphere metric/shell/SWE gates,
  the H&J background gate, the five pole gates (F0-F4) and the two global SWE
  forecast gates (F5/F6), and the three previously un-wired TFC scripts
  (`test_metric_reduction`, `test_freestream`, `test_terrain_stretched_smoke`).
  The four gates too slow for every-push CI (`test_sphere_gw` ~7.4 min,
  `test_hj_baroclinic` ~2.6 min, `test_hj_baroclinic_ridges` ~3.3 min,
  `test_hj_baroclinic_global` ~16.7 min) move to a new `workflow_dispatch`-gated
  `sphere-slow-gates` job. The two H&J global cases stay deliberately
  target-less (`diag = False`): they are production launch configs gated by
  physics thresholds, and mastering them would cost ~2.3 h per case. (sphere_global_golden_masters)
- CI runs the terrain test battery: transform units, h≡0 identity oracle,
  elliptic composition oracle, resting-atmosphere/uniform-flow/mountain-wave
  gates, the Agnesi smoke case, the Agnesi-vs-Smith analytic oracle, and the
  `test_agnesi_hydrostatic` golden master in the flow-solver suite. (terrain_phase9_ci_wiring)


Removed
^^^^^^^

- Removed dead code from the data-assimilation layer ahead of the ModelState
  repair: `localisation.py` (broken absolute import, zero callers), the legacy
  `utils.ensemble` wrapper class, `set_p2_nodes`/`set_rhoY_cells`,
  `params.converter`, and the `HSprojector_3t2D/2t3D` horizontal-slice
  projectors (the supported 2D x-y cases are native 2D since the ModelState
  refactor, so the projector guard could never fire). (da_dead_code)
- Removed the commented-out legacy `_update_solution_variables`/`_update_solution_arrays_jit`
  variant from `compute_advection.py` (the live function of the same name is unaffected). (dead_advection_variant)
- Removed the dead `utils/debug_helpers.py` module: its only symbol `pl_sol` had
  zero call sites anywhere in the tree, and the module pulled in a top-level
  `matplotlib.pyplot` import for nothing. (debug_helpers_dead_module)
- Removed dead code from the HDF5 writer: the never-called ``vortz`` / ``vorty``
  vorticity methods and ``dpress_dim`` (the last also referenced ``np.complex``,
  removed in NumPy >= 1.24), plus their commented-out dataset entries. (io_dead_vorticity_removed)
- Removed the unused `utils/operators/laplacian/lap2D_numba.py` (stencil-based
  2D laplacian): the production 2D elliptic path uses `lap2D_manual.lap2D_gather`
  exclusively (`@nb.njit`, via `implicit_euler._prepare_2d_system`), and nothing
  imported the module. The 3D path's `lap3D` kernel likewise remains numba-jitted. (lap2d_numba_dead_code)


Added
^^^^^

- Added the MWR-2022 OSSE reproduction tooling: `run_scripts/osse_mwr2022.py`
  drives the paper's five-run recipe per case (obs, truth, noda ensemble, LETKF
  with/without blending, plus ETPF variants) with fully seeded, regenerable
  observation files, and `run_scripts/osse_diagnostics.py` computes the
  ensemble-mean RMSE vs truth and ensemble spread per field over the
  assimilation window (CSV + PNG). The warm-bubble case gains the paper's seeded
  ensemble perturbation (`delth += 10*rand()`, truth seed 1234) — inert for the
  default deterministic run — and the travelling-vortex `logging.info` misuse
  in `sol_init` is fixed. (da_osse_scripts)
- The ETPF's optimal-transport step now uses POT (`ot.emd`) instead of the
  undeclared, unmaintained `pyemd` dependency, installable via
  `pip install "pybella[da]"` and import-guarded so the deterministic solver and
  the LETKF never need it. Both are exact solvers: on an ETPF-shaped random
  N=10 problem the transport plans agree to machine epsilon (max abs difference
  2.8e-17, identical objective cost), and `emd_with_flow`'s extra-mass penalty
  was irrelevant because both ETPF histograms sum to one. (da_pot_extra)
- Hughes & Jablonowski (2023) mountain-induced baroclinic wave (sphere
  follow-up), foundation: the Ullrich et al. (2014/2016) dry balanced base
  state as a reusable analytic module (`tests/ullrich_baroclinic.py`) — the
  closed-form temperature, pressure, density and gradient-wind zonal jet
  T/u/p/rho(phi, z) (Appendix B, dry variant q_v=0). pyBELLA being
  height-based, these evaluate directly with no pressure root-finding. The
  module also carries the two-ridge H&J topography (Eq. 1, Table 1) with
  analytic longitude/latitude gradients (required by `SphericalTerrainMap`).
  Gated by `test_hj_background.py`: surface T endpoints (310/240 K), the
  ~28 m/s midlatitude jet near 45 deg / 10 km, hydrostatic (dp/dz = -rho g)
  and gradient-wind self-consistency; the 2000 m ridge peaks, their 10%
  half-widths, longitude periodicity and finite-difference-checked gradients.
  Next: the field-mode hydrostatic reference + well-balanced sphere IC, the
  ridges wired through `SphericalTerrainMap`, and the multi-day run. (hj_ullrich_background)
- Reinstated nonhydrostatic<->hydrostatic dynamics-regime blending (thesis ch. 4.2)
  as a schedule-only mechanism, with no explicit conversion routine (matching the
  thesis-era RKLM_Python). The eos schedule (`physics/eos.py`) flips the blending
  parameter `alpha_w` (`is_nonhydrostatic`/`nonhydrostasy`) per step -- hydrostatic
  for the first `no_of_hy_initial` steps, then nonhydrostatic -- and `alpha_w` is
  wired into the two places that make the vertical-momentum equation degenerate to
  hydrostatic balance when it is zero: the elliptic vertical self-coupling
  `h22 = 1/nu_nh` (`numerics/coriolis.py` + JAX twin; drops the spurious `nonhydro`
  factor, giving the thesis coefficient of eq. 4.40/4.42) and the explicit
  vertical-momentum inertia (`numerics/explicit_euler.py`). Nonhydrostatic steps
  (`alpha_w = 1`) run the normal predictor->corrector; the hydrostatic step
  (`alpha_w = 0`) runs the thesis one-step recipe (sec. 4.2.3-4.2.4): two first-order
  predictor updates for an input-imbalance-independent balanced vertical momentum
  plus a second-order `pi` update, then the second-order corrector.
  Added the `test_blending_hydrostatic` regression case (imbalanced Skamarock-Klemp
  hydrostatic inertia-gravity wave) and wired it into CI. Reproduces the thesis
  sec. 6.2.3 w-probe suppression (centre-probe vertical-velocity rel-err ~0.022,
  matching the reported 0.024 within the harness normalisation). Bit-identical for
  `alpha_w = 1`. See `dev_notes/hydrostatic_blending.md`. (hydrostatic_blending)
- Added a physics oracle for the Baldauf & Brdar (2013) internal-gravity-wave
  case: `pybella/tests/baldauf_brdar_analytic.py` builds a numerically-exact
  solution of the linearised compressible Euler equations about the isothermal
  background (Bretherton-transformed constant-coefficient system, staggered
  vertical collocation, exact-in-time eigenpropagation per x-Fourier mode;
  energy drift ~1e-13) and evolves the simulation's own initial condition.
  `test_scripts/test_igw_analytic.py` gates the regression configuration
  against it (catching sign/dispersion/amplitude *wrongness*, complementing
  the golden masters which catch *change*) and writes ref/sim/diff PNGs.
  Refinement studies (dt 500->125 s, dx 20->10 km, f on/off) decompose the
  measured sim-vs-linear residual and surface a known solver limitation: in 2D
  runs the out-of-plane momentum receives only the implicit half of the
  Coriolis rotation (`explicit_euler.do_forward_step` skips the `rhow` row for
  `ndim == 2`), pinning that component's error at ~0.44 rel-L2. (igw_analytic_oracle)
- JAX backend, components 3+4 of 4 (advection + physics kernels; full-run
  validation): the per-sweep advection kernel (recovery + HLL upwind fluxes),
  the advective rhoY flux convolution, the Coriolis H^-1 apply/coefficients,
  and explicit diffusion now run on JAX under `ud.backend = "jax"`. The seam
  sits inside the numpy drivers (sweeps, flips, ghost fills and flux-difference
  updates are shared between backends), so flux-container semantics are
  bit-faithful. New env var `PYBELLA_BACKEND=jax` flips the backend for
  unmodified cases. Validation: driver-level equivalence on five states at
  magnitude-scaled 1e-12 (pure explicit arithmetic), and — end to end — the
  full regression suite runs on the JAX backend against the *same stored
  golden-master targets* as numpy (all fields ~1e-7 max-abs vs the 1e-5
  tolerance). Init-time code (hydrostatics, initial_pressure) intentionally
  stays numpy: it runs once and feeds both backends identically. (jax_backend_advection_physics)
- JAX backend, device-resident plan phase A: functional, jitted twins of every
  per-step boundary operation — ghost-cell fills (periodic/WALL no-gravity
  pads, the sequential two-layer hydrostatic gravity fill incl. terrain
  contravariant reflection and ATMOSPHERIC_EXTENSION), ghost-node fills (the
  periodic overlap exchange, reflect, quasi-2D degenerate broadcast), Rayleigh
  sponge damping (with func-mode forcing arrays), and wall-node rhs scaling —
  routed through the hybrid seams under `ud.backend = "jax"`. All gravity-fill
  index algebra, stratification-at-ghost-coordinates and metric slices are
  precomputed per orientation (canonical and sweep-flipped) at config build;
  the kernels are pure gathers + closed-form arithmetic. Validated per-BC on
  six states × flipped sweep orientations × sol-override ×
  compressible/incompressible at magnitude-scaled 1e-13, with np.pad-callback
  micro-tests (bit-faithful including the degenerate quasi-2D axis), plus the
  full 10-case regression gate on the hybrid backend. After this phase the
  hybrid step has no numpy kernels left — only orchestration. (jax_backend_boundary_fills)
- JAX backend, component 2 of 4 (elliptic solve): new `ud.backend` flag
  (default `"numpy"`) selects the solver for the semi-implicit pressure
  projection. With `backend = "jax"`, `implicit_euler.do_implicit_part` builds
  the JAX laplacian linop and solves with `jax.scipy.sparse.linalg.bicgstab`
  under scipy-matching convergence semantics (rtol 1e-5 / atol `ud.tol`);
  coefficient assembly, preconditioning and the rhs stay on the numpy path, so
  both backends solve bit-identical systems. Validated by residual
  certificates (the JAX solution satisfies scipy's stopping criterion under
  the numpy operator) and end-to-end `do_implicit_part` equivalence on six
  states at magnitude-scaled 1e-5 (the convergence-slack floor; the vortex
  case agrees at 3.6e-10). Default-path behaviour is unchanged. (jax_backend_elliptic_solve)
- JAX backend, component 1 of 4 (pure operators): `pybella.backends.jax_ops`
  mirrors `utils/operators` module-for-module (finite_difference, gradient,
  convolution, divergence, preconditioner, lap2D, lap3D) with functional,
  jit-clean twins of every kernel, validated against the numpy/numba
  implementation as a golden master at 1e-13..1e-14 (x64, magnitude-scaled)
  on realistic coefficient fields covering periodic/WALL/atmosphere
  boundaries and the terrain metric. jax is an optional dependency
  (`pip install pybella[jax]`); the canonical numpy path is untouched and the
  equivalence suite skips cleanly without jax. New CI job `jax-equivalence`
  runs the suite on Python 3.12 with an unpinned jax. (jax_backend_pure_operators)
- Device-resident JAX time loop (`ud.backend = "jax-device"`, env
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
  (code is hardware-agnostic). (jax_device_backend)
- Lifted the JAX boundary's field-mode + incompressible gravity-fill guard for
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
  is exact, and the case runs projection-free by default regardless. (jax_field_mode_incompressible_gravity_fill)
- Restored the Baldauf & Brdar (2013) linear internal gravity wave case from `archive/full_coriolis` as the regression test `test_igw_baldauf_brdar`. (restore_igw_baldauf_brdar)
- Restored a shallow-water regime demonstration as the regression case `test_swe_vortex`: a balanced SWE vortex (gamma = 2 gas-dynamics equivalence) ported from the legacy `balanced_shallow_water_2D` initial condition. (restore_swe_demo)
- Schär (2002) ridge golden-master case (`test_schaer_ridge`): two-scale
  orography (h0 250 m, a 5 km, lambda 4 km) under SLEVE coordinates, native 2D,
  with analytic gradients and the exact cos^2 smooth/residual split. New linear
  FFT mountain-wave oracle (`schaer_linear_analytic.py`, self-tested against
  Smith's closed-form drag) gates the wave field, drag and flux constancy, and
  the Gal-Chen-vs-SLEVE discriminator: spurious small-scale w aloft collapses
  ~17x under SLEVE (E_ss 0.004 vs 0.070) where the true lambda-scale response
  is evanescent-dead. Smoke + analytic tests wired into CI. (schaer_ridge_case)
- SLEVE vertical transform (Schär et al. 2002; Leuenberger et al. 2010
  exponent n): two-scale orography split `ud.orography_smooth` (large-scale)
  + residual against the total `ud.orography`, per-component sinh decay,
  analytic Jacobian — the first eta-dependent J through the operators.
  Selected via `ud.vertical_transform = SLEVETransform(s1, s2, n)`; the
  activation contract and the Gal-Chen path are untouched. Proven by
  transform identities, scale-separation property, and SLEVE-parametrized
  elliptic-composition + resting-atmosphere oracles (3D and native 2D). (sleve_transform)
- Spherical geometry, JAX twins (task 1a): the hybrid JAX backend now fills
  ghost cells for non-vertical-line (spherical) metrics.
  `backends/jax_ops/boundary.py` gains (1) the general free-slip WALL mirror
  (twin of `cell_boundary._mirror_momenta_general`: symmetric pad + per-side
  contravariant Cramer reflection with the local area normals), invoked at
  canonical orientation for the phi and thin-shell degenerate-r walls, and
  (2) the well-balanced radial gravity ghost fill (twin of the e_up branch of
  `_calculate_ghost_values`/`_assign_ghost_values`: hydrostatic rho/rhoY,
  tangential-velocity copy, and the beta·e_up momentum reconstruction that
  makes the rhoY wall flux exactly odd), in phys + sweep orientations. The
  non-vertical-line fast-fail guards are removed. Gated by
  `test_cells_sphere_*` in `test_jax_boundary_equivalence.py` and by
  `test_strang_sphere_{swe,gw}` (full sweep driver) in
  `test_jax_advection_equivalence.py` — JAX vs numpy at the ulp floor,
  including the mid-sweep (flipped) orientations. (sphere_jax_boundary)
- Spherical geometry, JAX twins (task 1b, part 1): the hybrid JAX backend
  now supports the sphere's Coriolis path. `backends/jax_ops/coriolis.py`
  gains a general-H⁻¹ kernel twin (`compute_coefficients_general` /
  `apply_inverse_general`) mirroring the numba Sherman-Morrison + (C11)
  alpha_w kernel for an arbitrary local up-direction `e`, and now consumes
  the spatially varying rotation field (`ud.coriolis_field`) and buoyancy
  rank-one term by reusing the numpy `role_components` / `_up_role_components`
  (cached, bit-identical). The `coriolis_field` fast-fail guard is removed.
  Gated by `test_coriolis_general_sphere_{swe,gw}` in
  `test_jax_advection_equivalence.py` (JAX vs numpy to 1e-13 on the TC2
  thin-shell and DCMIP-31 shell fixtures). (sphere_jax_coriolis)
- Spherical geometry, JAX twins (task 1, device-resident): the device-
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
  (TC2 device-vs-numpy at the Krylov floor). (sphere_jax_device)
- Spherical geometry, JAX twins: an in-process Williamson TC2 full-run gate
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
  the shared boundary-config no longer rejects them). (sphere_jax_fullrun_gate)
- Added the Straka density current (Straka et al. 1993) as regression case
  `test_straka` — the suite's first nonlinear, advection-dominated gravity case
  (256x32 at 200 m, dt = 4 s to t = 900 s; front position, peak winds and
  symmetry verified against the published benchmark). Includes a new explicit
  constant-coefficient diffusion module (`flow_solver/numerics/diffusion.py`,
  enabled per-case via `ud.diffusion` / `ud.diffusion_coeff`, off by default)
  implementing the benchmark's fixed K = 75 m^2/s on velocity and potential
  temperature. The case runs periodic in x: the x-WALL elliptic path is
  currently broken/untested (wall momentum zeroing in the nodal divergence
  handles the vertical axis only) — known limitation, to be fixed with the
  axial-agnosticity refactor. (straka_density_current)
- SWE <-> lake blending reinstated on the ModelState API (2D x-y layout):
  `do_swe_to_lake_conv` / `do_lake_to_swe_conv` ported, the orchestration's
  broken SWE branches repaired (undefined `flux`, unbound conversion flag,
  duplicate-`ud` latent TypeError on the continuous-blending call), and a new
  golden-master regression case `test_blending_swe` added (balanced SWE vortex
  with initial swe->lake->swe blend). (swe_lake_blending_reinstated)
- Stage F (poles) F0: pole-aware spherical metric build. ``SphericalShellMap``
  and ``SphericalTerrainMap`` gain a ``pole`` flag that folds ghost coordinates
  beyond |phi| = pi/2 through the pole (lambda -> lambda + pi, phi mirrored) so
  the ghost metric is a fold-copy of the far-side image cell (J > 0, no vector
  rotation); pole NODES carry J == 0 with ooJ := 0 and a stored ``pole_mask``.
  Adds ``BdryType.POLE`` with grid-init validation (even longitude cell count,
  2*pi longitude span, [-pi/2, pi/2] latitude). Gated by
  ``test_scripts/test_sphere_pole_metric.py``; channel cases bit-identical. (20260706_120000_sphere_poles_f0)
- Stage F (poles) F1: pole ghost exchange for cells and nodes. A ``BdryType.POLE``
  axis fills its ghost slabs by a pure index remap (longitude -> longitude + pi,
  latitude mirrored across the pole) with NO vector rotation — momenta are global
  Cartesian, so every field copies component-by-component. The source is always an
  interior cell (the far side of the pole), so the fill is order-independent within
  a sweep and works at any sweep-flipped orientation via the metric's tracked
  axes. Node pole rows are forced to their longitude-ring mean (single-valued at
  the pole). Gated by ``test_scripts/test_sphere_pole_exchange.py``. (20260706_130000_sphere_poles_f1)
- Stage F (poles) F2: pole-face flux closure. The lat-lon pole face has zero
  area in the continuum, but the discrete advective flux there is a symmetric
  O(dphi^2) residual that both pole-adjacent cells subtract with the same sign
  — a double-loss mass/tracer leak. The advection now zeroes the flux through
  the two pole faces (both the conservative HLL flux and the advecting mass
  flux), making the pole a no-flux edge; over-pole transport is carried by the
  longitude sweep. Gated by ``test_scripts/test_sphere_tc1_pole.py`` (Williamson
  TC1 cosine bell over the pole): J-weighted tracer conserved to machine
  precision, bell crosses the pole bounded, equatorial-rotation reference twin. (20260706_140000_sphere_poles_f2)
- Stage F (poles) F3: FFT-in-longitude polar filter. ``ud.polar_filter`` (a
  ``PolarFilter(phi_c, p)``) damps the CFL-violating zonal wavenumbers poleward
  of ``phi_c`` once per step with the classic transfer function
  r(k,phi) = min(1, [(cos phi/cos phi_c)/sin(pi k/N)]^p); it acts on the
  J-weighted conservative fields and leaves k=0 untouched, so every longitude
  ring's integral is conserved to roundoff (terrain-safe). The advective CFL is
  correspondingly relieved by cos(phi)/cos(phi_c) near the poles. ``None`` -> the
  solver is bit-identical. Gated by ``test_scripts/test_sphere_polar_filter.py``. (20260706_150000_sphere_poles_f3)
- Stage F (poles) F4: elliptic pole-ring collapse. Each (radius, hemisphere)
  lat-lon pole ring — one physical point spread over every longitude node —
  collapses to a single master pressure unknown via a Galerkin scatter/gather
  wrapped around the one-sided pole operator (``lap3D`` already treats the
  non-periodic latitude axis as a wall). The reduced system embeds in the
  full-size solve vector with non-master ring entries pinned to zero, so the
  BiCGSTAB plumbing is untouched; the rhs divergence zeroes the beyond-pole
  ghost momenta to match. Non-pole cases are bit-identical. Gated by
  ``test_scripts/test_sphere_pole_elliptic.py`` (resting global shell at rest,
  single-valued pole pressure, interior + ring-summed divergence to the floor). (20260706_160000_sphere_poles_f4)
- Stage F (poles) F5: Williamson TC2 on the full pole-to-pole sphere
  (``test_sphere_swe_tc2_global``, registered). The +-80 deg channel extended to
  |phi| <= pi/2 with POLE latitude walls, a pole-enabled shell map and the polar
  filter — the acid test that the pole ghost exchange, pole-face flux closure,
  polar filter and elliptic pole collapse compose in a real forecast. The steady
  zonal flow stays steady at the Krylov floor (l2(h) ~ 6e-6 over 10 steps),
  the resting global shell is exact, and the tangent-plane constraint holds to
  machine precision (surface_constraint is re-applied after the filter). Gated by
  ``test_scripts/test_sphere_swe_global.py``. (20260706_170000_sphere_poles_f5)
- Stage F (poles) F6: TC6 Rossby-Haurwitz on the full pole-to-pole sphere
  (``test_sphere_swe_tc6_global``, registered). Unlike TC2, TC6 carries genuine
  longitude structure at all latitudes, so it exercises the pole ghost exchange
  and polar filter under real dynamics. Over a short run the wave stays bounded,
  J-weighted mass is conserved to machine precision (pole-face closure under
  dynamics), and the RH-4 zonal-wavenumber signature is preserved. The 7-day
  phase-speed validation runs from run_scripts. Gated by
  ``test_scripts/test_sphere_swe_tc6_global.py``. (20260706_180000_sphere_poles_f6)
- Stage F (poles) F7a: JAX HYBRID (``backend='jax'``) twins of the pole
  machinery. The lat-lon pole fold is a pure functional index remap of the
  phi-ghost slabs (cells + nodes, all 6 fields, no vector rotation), reusing the
  same ``common.pole_source_indices`` maps as numpy; the node fill adds the
  pole-row longitude ring-mean. The elliptic pole-ring collapse is twinned as a
  Galerkin scatter/gather wrapped around the JAX ``lap3D`` matvec (one-sided pole
  rows for free — the kernel already slab-zeroes the non-periodic phi axis); the
  ring is zeroed by a mask multiply (bicgstab's ``custom_linear_solve``
  double-transposes the operator, which an integer scatter-SET cannot survive).
  The polar filter runs unchanged on the shared ``time_update`` path. Gates:
  pole cell/node fills numpy-vs-jax to machine precision (canonical + phi-sweep),
  the collapsed ``do_implicit_part`` on a global shell to machine precision, and
  ``test_sphere_tc2_global_{stepper_bit_identical,hybrid_reproduces_numpy}`` —
  the pole-to-pole TC2 stepper reproduces numpy bitwise (projection off) and at
  the initial-projection Krylov floor (projection on). ``jax-device`` still
  fast-fails on POLE (device twins land in F7b). (20260706_190000_sphere_poles_f7a)
- Stage F (poles) F7b: JAX DEVICE-RESIDENT (``backend='jax-device'``) twins of
  the pole machinery. The pole ghost fold is wired into the device ``_ghost_fill``
  (the shared ``BoundaryConfig`` pole_cell_fill, canonical + phi-sweep); the
  elliptic pole-ring collapse is applied inside the device elliptic closure
  (Galerkin scatter/gather around ``jax_lap3D._lap3D``, rhs gathered, master
  scattered back post-solve); and — the biggest device-specific piece — the FFT
  polar filter runs INSIDE the jitted step loop (``device.run_window`` returns
  before the numpy filter attach point), followed by the tangent-plane surface
  constraint and a ghost refill, mirroring ``time_update.do``. The longitude-CFL
  cap is applied in the device dt control. All are strictly guarded no-ops when
  there is no POLE axis / no ``polar_filter``, so every existing device golden
  master (TV, Straka, 3D vortex, Agnesi, channel TC2/SWE) stays bit-identical
  (the dispatch-guard proof). Gates: ``test_window_sphere_swe_global``
  (device-vs-hybrid ~1e-13, compile-count 2) and
  ``test_sphere_tc2_global_device_reproduces_numpy``. The device pole path ran an
  adversarial ``pybella-jax-reviewer`` pass over the fill + collapse + filter
  twins. (20260706_200000_sphere_poles_f7b)
- Stage F (poles) F8: Hughes & Jablonowski (2023) on the FULL pole-to-pole
  sphere — the +-80 deg channel cases extended to |phi| <= pi/2. Two registered
  cases ``test_hj_baroclinic_global`` (flat Ullrich background) and
  ``test_hj_baroclinic_ridges_global`` (the two midlatitude ridges), each
  subclassing its channel case with a ``pole=True`` map (deep ``SphericalShellMap``
  / ``SphericalTerrainMap``), ``BdryType.POLE`` latitude walls and the FFT polar
  filter (phi_c = 70 deg, poleward of the 45 deg N ridges). Two channel
  workarounds are lifted in the global variants (kept in the channel via a
  ``ud.phi_clip`` default): the 89.5 deg Ullrich cos^K clamp (a node now sits
  exactly at the pole, where the jet vanishes and the pressure is finite) and the
  "phi grid >= 32" constraint (ghosts fold instead of overshooting the pole).
  Gate ``test_scripts/test_hj_baroclinic_global.py``: the flat background stays
  steady (jet ~27.4 m/s, meridional adjustment 0.13 m/s — DROPS from the channel's
  0.18, no +-80 deg walls); the ridges initiate a longitude-localised meridional
  response (~1.57 m/s at 140.6 deg E, ~12x the flat background, jet bounded) with
  J > 0 at the ridge latitudes (the /h_ref nondimensionalisation holds); and the
  device-resident backend reproduces numpy at the Krylov floor. The multi-day
  production run stays GPU-gated (pt-3 protocol). Stage F (poles) is COMPLETE. (20260706_210000_sphere_poles_f8)
- Native-2D Agnesi equivalence proof: the 2D terrain path reproduces the
  quasi-2D 3D golden-master run to the pre-existing lap2D/lap3D wall-convention
  floor (~10% on the small wave fields, solver-tolerance independent, same gap
  as a flat wall-bounded impulse) and — the strong gate — passes the Smith
  (1980) analytic oracle with the same thresholds and near-identical
  calibration (w 0.397, u' 0.432, drag 0.981, flux 3.9%) at ~5x less runtime. (agnesi_2d_equivalence)
- Added the axis-geometry module `pybella/utils/axes.py` (role-space convention:
  cyclic permutation mapping (h1, v, h2) roles onto array axes, vertical-axis
  accessors, slab/profile/permutation helpers) with unit tests, plus the
  bit-for-bit H5 run comparator `test_scripts/compare_h5_runs.py` used to gate
  the pure-refactor phases of the axial-agnosticity work. No behaviour change. (axial_phase0_axes_module)
- Axial-agnosticity Phase 5: the permutation oracle
  (`test_scripts/test_permutation_oracle.py`, CI-wired) proves the endgame —
  the 2D internal-long-wave reference (gravity, stratification, walls, full
  Coriolis) embedded as z-vertical (gravity_direction=2) and x-vertical
  (gravity_direction=0) quasi-2D 3D twins reproduces the sigma-mapped 2D
  fields through full solver steps to <=1e-6, with exact uniformity along
  the degenerate axis and the Coriolis pseudovector mapping produced
  automatically by the role-based configuration. Also fixes a latent bug the
  oracle exposed: SpaceDiscr stored `dxyz/ig/ic/stride` as class-level shared
  arrays, so two grids coexisting in one process corrupted each other; they
  are now per-instance. Bit-identical for all existing cases. (axial_phase5_permutation_oracle)
- Hughes & Jablonowski (2023) mountain baroclinic wave, pt 3 (device path):
  a jax-device correctness gate for the ridge case's production path
  (`test_scripts/test_hj_baroclinic_ridges.py::test_ridges_device_reproduces_numpy`,
  jax-skip guarded). The ridge case is the UNION of two already-device-validated
  paths — the general non-vertical-line SPHERE metric (general e_up buoyancy,
  general H^-1, the constant `coriolis_field`, free-slip phi walls; TC2) and
  TERRAIN-following coordinates (agnesi) — now exercised together with radial
  GRAVITY, the terrain tilt, and a field-mode COMPRESSIBLE HydroState. Being
  compressible with no initial projection, it never hits the JAX boundary's
  field-mode+incompressible guard, so the numpy-only projection path is
  untouched. The device backend reproduces numpy at the per-step bicgstab
  Krylov / ulp floor (3 steps, 32x12x32: rho/rhoY/rhoX ~3e-7 abs, momenta/rho
  ~7e-5, p2_nodes ~3e-6 relative since p2 ~ O(1/Msq) ~ 700, all finite,
  device compile count 2) — and ran ~9x faster than numpy even on CPU. The
  multi-day full-planet production run itself is H100-class (no GPU on the dev
  box; sphere JAX was H100-validated) and is documented as a runnable protocol
  in dev_notes: `PYBELLA_BACKEND=jax-device pybella -ic test_hj_baroclinic_ridges
  -N 1` (target-less, `diag=False`) at production resolution, acceptance =
  qualitative Rossby wave train vs the paper's dry FV/SE panels + self-convergence. (hj_pt3_device_validation)
- Hughes & Jablonowski (2023) mountain baroclinic wave, pt 2: the two
  midlatitude ridges (the wave TRIGGER) via a terrain-following spherical map
  plus the well-balanced "adjusted" background (`tests/test_hj_baroclinic_ridges.py`,
  smoke gate `test_scripts/test_hj_baroclinic_ridges.py`). The ridges (Eq. 1,
  h0 = 2000 m at 72 E / 140 E, 45 N) enter as GEOMETRY: a
  `SphericalTerrainMap` (Gal-Chen radial coordinate, orography nondimensionalised
  to the map's `z / h_ref` units) with no SWE bottom-topography source terms.
  The balance is the surface-pressure adjustment (Eq. 2), which is exactly the
  Ullrich pressure profile (Eq. B4) sampled at the surface height
  `p_s(lambda, phi) = pressure(phi, z_s)`; more generally the whole adjusted
  state is the analytic base state sampled at the terrain-following height
  `z = r - a = Z(eta, h(lambda, phi))`. pt 1's `sol_init` is written entirely in
  `metric.height`, so it reproduces the adjusted state unchanged once the map
  carries the ridges — the case is pt 1 with the map swapped and the radial
  axis carrying eta in [0, depth] (verified: ridge-tip surface pressure 779 hPa
  vs the paper's ~773 hPa; J = r^2 cos(phi) stays positive at the ridge
  latitudes only once the orography is nondimensionalised). The adjusted state
  is well-balanced but not PERFECTLY so, and the residual near the ridges is
  the intended trigger. Gate: over a coarse ~40 min run the ridges drive a
  meridional-wind response ~8x the flat-background pt-1 adjustment (1.45 vs
  0.18 m/s) LOCALISED at the ridge centres (peaks at 73 / 141 E), with the jet
  bounded and nothing blowing up — the decisive ridge-triggered-initiation
  signature. Next: the multi-day full-planet device run (pt 3). (hj_pt2_ridges)
- Hughes & Jablonowski (2023) mountain baroclinic wave, pt 1: the well-balanced
  Ullrich base state on the deep spherical shell as a pyBELLA initial condition
  (no topography) with a steadiness gate (`tests/test_hj_baroclinic.py`,
  `test_scripts/test_hj_baroclinic.py`). The 3D balance is split the way the
  discretisation carries it: the VERTICAL hydrostatic balance goes into a
  z-only equatorial-Ullrich field-mode HydroState reference (built by the
  terrain-quadrature branch of `hydrostatics.integrated_state` from a custom
  `ud.stratification`), while the MERIDIONAL pressure structure goes into the
  Exner perturbation `p2_nodes` — the momentum pressure-gradient force uses only
  grad(p2), so the meridional gradient the Coriolis force on the jet must
  balance has to live there (the 3D analogue of TC2's geostrophic depth in p2).
  Full Coriolis is the constant embedded rotation vector `2*Omega_nd*(0,0,+1)`
  (the mirrored-embedding pseudovector flip). The registered case
  `test_hj_baroclinic` runs `SphericalShellMap` (full-size a, deep ~30 km shell,
  lambda-periodic, radial gravity, +-80 deg latitude free-slip walls). Gate:
  with no ridges the background stays steady — the ~28 m/s midlatitude jet holds
  to <0.05%, the geostrophic adjustment saturates the meridional wind at
  ~0.18 m/s (0.7% of the jet) and then plateaus, nothing blows up; a wrong
  Coriolis sign/factor would drive O(jet) meridional wind within a few steps.
  Next: the ridges via `SphericalTerrainMap` + the Eq. 2 adjusted balance, then
  the multi-day full-planet device run. (hj_pt1_wellbalanced_ic)
- NEDAS Phase D2: ensemble-vmapped GPU forecasts. `device_batch.run_window_batch`
  advances the whole ensemble in lockstep on device — batch-min host-controlled
  dt, the full compiled step (bicgstab included) vmapped over the member axis,
  optional member-axis sharding across GPUs (`ens_devices`) — wired into the
  `ens_run_strategy: batch` branch of `PyBellaModel.run()` for JAX device
  backends (`ens_batch_mode: vmap|loop`; `loop` is the same-dt gate comparator).
  Gate: `test_scripts/test_nedas_vmap_gate.py`. Design + validation:
  dev_notes/nedas_interface.md Phase D. (nedas-interface-d2)
- NEDAS Phase D4: blending + CFLfixed on the jax-device backend. Blend windows
  segment at the conversion steps: the device drivers run the same
  `schemes.prepare_blending` as the numpy loop (host round-trip only on
  conversion steps), each regime compiles its own static step variant
  (`is_compressible` joins the step-cache key), and the psinc variant exports
  the predictor half-time pressure the psinc→comp trial extraction reads. The
  pinned `initial_blending: True` configs, EnDAB, and the bubble `CFLfixed`
  case now run on GPU. Fixed (both sides): the blend trial integration is
  capped at exactly one step — previously a last-ULP dt tie could append a
  spurious dt≈0 step whose degenerate half-time pressure became the whole
  blended dp2n (backend/BLAS-dependent per-member results). (nedas-interface-d4)
- NEDAS Phase E: the first true 3D DA — igw3d OSSE through NEDAS. The adapter
  gains a 3D branch (horizontal (x,z) grid + vertical levels with per-level
  read/write and interface z_coords), PyBellaObs a seed-778 volumetric obs
  convention with a VarCov err-std floor, the internal-long-wave IC a 3D
  branch + seeded theta'-wave member perturbations (2D path bit-identical;
  oracle-gated), and the OSSE tooling igw3d generation, ndim-general
  diagnostics with the w=rhov/rho blending probe, and dev/prod configs.
  Self-anchored gates (no native 3D pipeline exists):
  test_igw3d_ic_oracle.py, test_nedas_obs3d_selfcheck.py, and the igw3d
  variant of test_nedas_vmap_gate.py. Findings + ladder results in
  dev_notes/nedas_interface.md Phase E. (nedas-interface-e1)
- Phase N0 of the NEDAS interface: skeleton `pybella.interfaces.nedas` adapter
  (PyBellaModel/PyBellaObs stubs + registry hook), `pybella[nedas]` optional
  extra pinned to NEDAS==1.2.0, and a travelling-vortex OSSE run config
  (`run_scripts/nedas_tv_osse.yml`). Recon findings + adapter design in
  dev_notes/nedas_interface.md. (nedas-interface-n0)
- Phase N1 of the NEDAS interface: working `PyBellaModel`/`PyBellaObs` adapter
  (in-memory per-member ModelState forecasts, native-parity observations),
  `run_scripts/nedas_run.py` launcher + `nedas_osse_diagnostics.py` exporter,
  and `test_scripts/test_nedas_obs_parity.py` (obs byte-identical to the frozen
  native pipeline). TV noda/EnDA/EnDAB OSSEs run end-to-end through NEDAS. (nedas-interface-n1)
- NEDAS N3 kickoff: rising-bubble OSSE configs, TopazDEnKF example config,
  JAX-backend support verified through NEDAS, and a bit-inert
  `ens_run_strategy: batch` hook in PyBellaModel (the future vmap insertion
  point for ensemble-parallel JAX forecasts). (nedas-interface-n3)
- NEDAS adapter: p2_nodes joins the DA state (cell-centred view + increment-only
  node write-back; native cell/node analysis split via impact_on_state) and
  PyBellaObs generalises to node-grid observations. Fig-10 all-quantities TV
  OSSE configs added; obs byte-parity test extended to the 5-attribute set. (nedas-interface-p2)
- Spherical geometry, stage B: spatially varying Coriolis —
  `ud.coriolis_field` (callable of the Cartesian coordinates returning the
  3 Cartesian rotation components, evaluated once per grid and cached)
  feeds the H^-1 kernels and the explicit forward step as per-cell fields;
  the legacy scalar `ud.coriolis_strength` path is bit-identical.
  `SphericalShellMap.traditional_coriolis` provides f(phi) e_r; the JAX
  backends fast-fail on rotation fields until the sphere SWE JAX stage. (sphere_pt2_coriolis_field)
- Spherical geometry, stage A2: general curvilinear free-slip walls — the
  ghost-cell fill mirrors the CONTRAVARIANT momentum triple with the local
  area normals (per-cell 3x3 solve; exact wall-flux cancellation) for
  non-vertical-line metrics, in both the no-gravity WALL path and the
  gravity path (tangential-velocity copy + `e_up` normal reflection,
  `h_v`-based ghost spacing). Vertical-line maps keep the legacy slope-term
  recipe verbatim (bit-identity); the JAX ghost fill fast-fails on
  spherical metrics until the SWE stage. (sphere_pt1_general_walls)
- Spherical geometry, stage D5: DCMIP-31-style nonhydrostatic gravity wave
  on the small-planet compressible shell (`test_sphere_gw`) — a
  latitude-independent, longitude-periodic potential-temperature
  perturbation on the isothermal (Baldauf-Brdar) background whose
  large-radius equatorial slice reduces to the planar B&B channel. Fixes
  the general (curved-metric) free-slip gravity wall to be
  FREESTREAM/WELL-BALANCED: the wall ghost now reflects the rhoY flux
  `Y*(N_v.m)` as exactly odd (using the image cell's `N_v.e_up`, not the
  source cell's), so the wall mass flux cancels to roundoff and J-weighted
  mass/energy are conserved to machine precision under dynamics — where the
  prior up-velocity reflection left an O(dz^2) leak. Confined to the
  `vertical_line=False` path; vertical-line terrain and uniform-Cartesian
  walls stay bit-identical (all terrain/shell/SWE golden masters green).
  Stage D physics gates complete: J-weighted conservation machine-exact,
  2nd-order self-convergence (Richardson 3.92), and the large-radius limit
  converging to the planar Baldauf-Brdar linear oracle (zonal-field rel-L2
  0.187 -> 0.098 -> 0.051 as X = 125 -> 62.5 -> 31.25). (sphere_pt6_gravity_wave)
- Spherical geometry, stage D: radial gravity on the true-radius 3D shell —
  buoyancy and the alpha_w (hydrostasy) structure act along the local up
  `metric.e_up` (explicit forward step, implicit explicit-part, rhoX
  stratification coupling), the H^-1 Coriolis/buoyancy matrix generalizes
  to an arbitrary up direction (Sherman-Morrison closed form + the thesis
  (C11) alpha_w prefactor `G = H^-1 - (1-a)(I-ee^T)H^-1(ee^T)`, reducing
  to the legacy kernel term-by-term <= 1e-14), and the isothermal
  hydrostatic state integrates along radial columns (height/h_v). Gate:
  the resting shell stays at rest to the Krylov floor (|u| ~ 2e-12),
  including an a = 5-scale-heights small planet. Legacy Cartesian cases
  keep the (C11) kernel verbatim (bit-identity); SWE-channel golden
  masters regenerated (the general-kernel dispatch moves their initial
  projection within its solver-floor equivalence class). (sphere_pt4_radial_gravity)
- Spherical lat-lon geometry, stage A: `SphericalShellMap` (Tier-3
  `CurvilinearMap`; lambda-r-phi grid, Cartesian momenta, metric-as-data),
  `ud.curvilinear_map` activation in `grid_init`, `MetricFields`
  generalized-altitude data (`height`, `e_up`, `h_v`, `vertical_line`),
  and the `test_sphere_metric.py` oracle suite (duality, metric identity,
  frozen-shell radial-defect contract, gradient/divergence oracles,
  elliptic composition identity on the spherical channel). (sphere_pt0_shell_map)
- Spherical geometry, stage C: shallow water on the spherical channel —
  Williamson TC2 (steady geostrophic zonal flow) and TC6 (Rossby-Haurwitz
  R=4) as thin-shell gamma=2 SWE regression cases on the lat-lon channel
  (walls at +-80 deg), with the tangent-plane momentum constraint
  (`ud.constrain_to_surface`) and the traditional-approximation Coriolis
  field. Pins the mirrored-embedding rotation-vector convention
  (`rotation_axis_cart`): TC2 discrete balance holds at ~5e-8 after 10
  steps (wrong sign: ~5e-4), the resting shell is exact, and tangency
  stays at machine zero. (sphere_pt3_swe_channel)
- Spherical geometry, stage E: terrain-following coordinates on the sphere
  — `SphericalTerrainMap` composes the existing `VerticalTransform`
  (Gal-Chen) into the radial coordinate, r = a + Z(eta, h(lambda, phi)),
  with lambda-wrapped orography, radial coordinate lines (t_eta || e_r)
  and the new `up_direction` map hook keeping GRAVITY radial while the
  coordinate surfaces tilt with the slope. Gates: h == 0 reduces to the
  shell map; duality + 2nd-order metric identity over a wavy hill (the
  freestream tripwire); the resting isothermal atmosphere over a spherical
  Gaussian-belt hill stays at rest to the solve floor. (sphere_pt5_terrain)
- New golden-master regression case `test_agnesi_3d`: an isolated circular 3D
  Agnesi bell (h0 = 100 m, a = 10 km, N = 0.01 1/s, U = 10 m/s, Gal-Chen,
  Rayleigh sponge above 12 km) on a 64x32x64 grid, 10 spin-up steps. First
  case where the orography depends on BOTH horizontal coordinates, locking the
  G2 != 0 terrain dynamics (second-slope contravariant fluxes, pressure-map
  column, elliptic cross terms, bottom BC) as a golden master — closes
  terrain-following limitation #4. Not wired into CI (full-3D elliptic solves,
  ~5 min/run); run locally via `pybella -ic test_agnesi_3d -N 1`. (terrain_agnesi_3d_golden_master)
- Agnesi hydrostatic mountain-wave case (`test_agnesi_hydrostatic`): N=0.01 1/s,
  U=10 m/s over a 100 m witch-of-Agnesi hill with a=10 km in Gal-Chen
  terrain-following coordinates, RAYLEIGH sponge above 12 km. Golden-master
  regression run (15 steps) plus a Smith-1980 analytic oracle
  (`test_agnesi_analytic`): quasi-steady wave field vs the linear hydrostatic
  solution — momentum flux matches the analytic wave drag to 2%, constant with
  height to 4% below the sponge. (terrain_phase8_agnesi_case_oracle)
- Terrain-following coordinates, phase 0: `VerticalTransform`/`GalChenTransform` +
  `MetricFields` scaffolding (`flow_solver/discretisation/terrain.py`), built in
  `grid_init` and attached as `elem.metric`/`node.metric` (`None` without
  `ud.orography` — uniform-Cartesian path untouched). Transform unit tests, h≡0
  identity oracle, and a quasi-2D mountain-wave smoke case (`smoke_agnesi`)
  through the full-tensor 3D elliptic path. (terrain_phase0_transform_scaffolding)
- Generalized terrain metric (tfc pt 1): `MetricFields` now carries the Klein
  area normals `N` (outer index = array axis, rotating with sweep flips; inner
  index = fixed Cartesian component) and physical coordinates `x`, synthesized
  bit-exactly from the vertical-line map when not supplied. New `CurvilinearMap`
  base + `build_metric_fields_from_map` general-path builder (Tier-2 capable,
  J > 0 guarded, effective slopes for legacy consumers); device config mirrors
  the normals. (tfc_phase1_general_metric_fields)
- Phase-0 reduction contract for the curvilinear (Klein) metric generalization:
  `test_scripts/test_metric_reduction.py` pins the vertical-line tangents/normals,
  the duality `N_i . t_j = J delta_ij`, `det N = J^2`, the bit-exact flux-reduction
  contract against `_metric_contravariant_fluxes_jit` (3D and 2D), and a Tier-2
  stretched-map tangent fixture with the effective-slope identity. (tfc_phase0_metric_reduction_contract)


Changed
^^^^^^^

- Routed backend detection through the existing ``backends.is_jax_backend``
  helper. The literal ``getattr(ud, "backend", "numpy") in ("jax", "jax-device")``
  check was copy-pasted inline across eight hot/parity-path modules (cell/node/
  common/rayleigh boundaries, advective flux, diffusion, coriolis, implicit
  elliptic); they now call ``is_jax_backend(ud)`` — the single place the
  backend-name set is defined. Behaviour-identical (the helper *is* that
  expression): numpy bit-identical, jax/jax-device green. (backend_dispatch_helper)
- Flattened the blending orchestration: named predicate helpers
  (`_window_start_conversion_due`, `_full_blend_due`, `_initial_blend_phase`)
  replace the nested conditionals, the field-order-dependent `ModelState`
  iterator unpack is gone, and the dead `debug` threading is dropped.
  Bit-identical on all three blending masters and the fast set. (blending_orchestration_flatten)
- Hardened the `test_blending_warm_bubble` regression case: the IC now carries a
  large Exner-pressure blob (amplitude 1.0, ~10x the slow signal, momenta
  untouched) that the initial comp->psinc blend must absorb in one step, and the
  historical loose tolerances (1.0 on momenta, 1e-4 on p2 increments) are
  replaced by the 1e-5 defaults on all 7 fields. Rationale: the old
  weakly-imbalanced IC plus loose tolerances could not detect the
  discarded-conversion defect (momenta shift ~0.1, final p2 increment diff
  ~1e-5) — verified by running the pre-fix code against the new master, which
  now fails at rho max-abs 1.1e-3 vs 1e-5 while the fixed code reproduces to
  ~1e-7. A uniform pi offset was shown to be gauge-inert (pi acts only through
  gradients), so the blob is the minimal genuine imbalance; blob absorption and
  the blended trajectory are validated run-vs-run against the frozen paper-era
  code (archive/localdab, cross-code agreement 4e-4). Golden master regenerated
  deliberately; blending-swe/hydrostatic and the fast set are untouched. (blending_warm_bubble_blob_ic)
- Factored the remaining copy-pasted regression-case idioms into
  `tests/case_setup.py`: `do_initial_projection` (the ~35-line incompressible
  initial-projection block, duplicated verbatim across the travelling-vortex, SWE
  and 3D-Coriolis cases) and `mirror_centers` (the periodic nearest-image vortex
  centre, duplicated 3×). The unused `T_from_p_rho` helper (defined in 5 cases,
  called in none) and the local `class obj` attribute-bag (replaced by
  `types.SimpleNamespace` inside the helper) were removed, along with the imports
  they orphaned. Test-only; FAST/affected cases bit-identical at tol 0. (case_setup_initial_projection)
- Ported the DA layer's member access onto the ModelState API: new
  `data_assimilation/ensemble_access.py` is the single place that knows on which
  container a DA attribute lives (`CellSolField` vs `NodePressureField`),
  replacing the pre-refactor `results[:, loc, ...]` container-index convention
  throughout `analysis.py`, `letkf.py`, `etpf.py` and `utils.py`. The `dap.loc`
  index dict and magic ghost-pad widths are gone (`elem.igx/igy` now). The
  LETKF (Hunt et al. 2007) and ETPF algorithmic cores are verified unchanged
  against the pre-refactor reference (`7f0b676~1`) modulo black formatting. (da_ensemble_access)
- Split the jax-device host window driver (`run_window`, CFL/dt host control,
  support guard, forcing eval, compile cache) out of `device_kernels.py` into
  `backends/jax_ops/device_loop.py`; `device_step` re-exports are unchanged. (device_loop_split)
- Split the 1273-line ``backends/jax_ops/device_step.py`` into ``device_config``
  (the host-side ``DeviceConfig`` + laplacian gather/mask-plan builders),
  ``device_state`` (host↔device 7-leaf-pytree marshalling), and ``device_kernels``
  (the traced substeps, the compiled ``make_step``, and the host ``run_window``
  loop). ``device_step`` keeps its design docstring and becomes a re-export shim,
  so ``device_step.run_window`` is unchanged. No arithmetic moved: jax-device
  output is bit-identical to the pre-split run and the inline gate stays green. (device_step_split)
- Routed the last six regression cases (agnesi, Schär, IGW B&B, blending warm
  bubble, unstable Lamb, 3D-Coriolis travelling vortex) through
  `tests/case_setup.make_diag_state` instead of constructing `DiagnosticState`
  inline, so all eleven cases now share the one helper that centralises the
  `Nx=inx-1 / Ny=iny-1 / steps=[stepmax-1]` offsets. Case-specific tolerances,
  `time_increment`, the f-string `test_name` and the rationale comments are
  preserved verbatim as forwarded keywords. Proven value-identical: `vars(diag_state)`
  byte-identical for all six, plus bit-identical solver output at tol 0. (diag_state_helper_migration)
- Reinstated the ETPF rejuvenation term (`+ delta * randn(...)` after the
  transport step), which was commented out in the reference code and caused
  severe ensemble collapse (spread ~5e-4 at rejuvenation_factor 0.001) in the
  MWR-2022 OSSE reproduction. (etpf_rejuvenation)
- Named the terrain hydrostate quadrature grid constants in `physics/hydrostatics.py`
  (`_HYDRO_QUAD_MIN_POINTS = 2048`, `_HYDRO_QUAD_POINTS_PER_CELL = 16`) with a comment,
  replacing the bare `max(2048, 16 * int(elem.sc[vv]))` magic numbers. Values
  unchanged; bit-identical at tol 0. (hydrostatics_quad_constants)
- Split the 800-line monolithic ``utils/io.py`` into an ``utils/io`` package with
  cohesive submodules — ``writer`` (the ``hdf5`` writer + ``initialise`` bootstrap),
  ``restart`` (``read_input`` / ``sim_restart`` / ``fn_gen``), ``debug`` (the debug
  writers), and ``cli`` (``get_args`` / ``init_logger`` / ``mkdir_p``). The public
  names are re-exported from the package ``__init__`` so every ``from ..utils import
  io`` call site is unchanged; output is bit-identical. (io_package_split)
- Replaced the bare ``21.69`` literal in the time-stepper's ``CFLfixed`` override
  with a documented module constant ``_CFLFIXED_DT_SECONDS`` (the legacy warm-bubble
  fixed timestep in seconds). Value unchanged. (prestep_cflfixed_constant)
- Split the 793-line `interfaces/dynamics_blending/schemes.py` god-module into a
  `schemes/` package: `blending` (the `Blend` interface), `comp_psinc`, `swe_lake`,
  `hydro_nonhydro` (commented hydro block kept verbatim as the reinstatement
  reference) and `orchestration` (the per-timestep blending calls). Public names are
  re-exported from `schemes/__init__.py`, so `time_update.py` and `prepare.py` call
  sites are unchanged. Mechanical move only — FULL repro gate bit-identical at tol 0
  (all 11 golden masters, incl. the blending warm bubble). (schemes_package_split)
- Factored the two repeated, gotcha-prone ``UserData`` idioms in the regression
  cases into ``tests/case_setup.py``: ``build_bdry`` (the per-instance,
  object-dtype boundary-type triple — never a shared class attribute) and
  ``make_diag_state`` (centralises the ``Nx = inx - 1`` / ``Ny = iny - 1`` /
  ``steps = [stepmax - 1]`` offsets while forwarding case-specific keywords).
  Behaviour-preserving: ``vars(UserData())`` is byte-identical for all 11 cases
  and the regression gate stays bit-identical. (test_case_setup_helpers)
- Restructured outermost looping in main() (23a7817db225fe2188bcf3045cfdc863435b6888)
- Axial-agnosticity Phase 3: hydrostatic state, vertical profiles, and
  boundary handling are axis-generic — `States` profiles carry their vertical
  axis (`expand_profile`), `integrated_state`/`analytical_state` integrate
  along the configured vertical, the nodal-divergence wall zeroing loops over
  all WALL/RAYLEIGH axes (fixing the previously broken x-WALL elliptic path),
  the gravity ghost-cell handler threads the physical vertical axis (also
  fixing a latent 3D bug where it read `gravity_strength[2] = 0` during
  advection sweeps), `_set_boundary` mirrors the wall-normal momentum of the
  actual wall axis, and the quasi-2D nodal broadcast generalises to any
  degenerate axis. Bit-identical for all existing cases (verified at tol=0). (axial_phase3_boundaries_hydrostates)
- Axial-agnosticity Phase 1: `gravity_direction` is now a configuration input
  (default 1 = y-vertical, validated; 2D runs require 1), `gravity_strength` and
  `coriolis_strength` are computed role-based from it, and all seven hardcoded
  `ud.gravity_strength[1]` reads route through `axes.vertical_axis(ud)`.
  Bit-identical for all existing cases (verified at tol=0 on full run outputs). (axial_phase1_gravity_config)
- Axial-agnosticity Phase 4: the 3D elliptic operator now carries the full
  H^-1 tensor coefficients (C_ij = (Gamma^-1 P Theta) h[role(i),role(j)], the
  same H^-1 the momentum correction applies), replacing the legacy ad-hoc x-z
  `corrf` cross terms; the 3D preconditioner uses the C_ii diagonals. With
  identity H^-1 the operator reduces bit-exactly to the previous one (Oracle
  A), and with full Coriolis a y-uniform 3D solve now matches the trusted 2D
  path under the pseudovector axis mapping to ~1e-6 (new
  test_scripts/test_3d_coriolis_oracle.py). That oracle also exposed and
  fixed the implicit-side sibling of the 2D out-of-plane Coriolis defect (the
  pressure-correction w-row was ndim==3-only). Regenerated golden masters:
  travelling_vortex_3d_coriolis (operator upgrade), igw_baldauf_brdar,
  internal_long_wave, unstable_lamb (w-row fix; the igw analytic-oracle
  out-of-plane error improves 0.060 -> 0.045). All other cases bit-identical. (axial_phase4_lap3d_full_tensor)
- Axial-agnosticity Phase 2: Coriolis and explicit dynamics are now
  role-canonical — `multiply_inverse_terms` binds (wh1, wv, wh2) and the
  momentum components in (h1, vertical, h2) role order via the cyclic axis
  permutation (njit kernels unchanged); new `compute_inverse_coefficients`
  exposes the cached H^-1 fields for the upcoming tensor elliptic operator;
  buoyancy and the rhoX stratification coupling act on the configured vertical
  momentum; the explicit momentum rows are written in role symbols with
  expression trees preserved; the advection pwchi special case keys on the
  vertical axis. Bit-identical for all existing cases (verified at tol=0). (axial_phase2_role_canonical)
- Axial-agnosticity Phase 6 wrap-up: Straka now runs its faithful free-slip
  wall configuration on all boundaries, exercising the repaired x-WALL
  elliptic path (target regenerated for the wall config); a z-vertical
  plumbing smoke case (`smoke_zvert` + `test_scripts/test_zvert_smoke.py`)
  pushes `gravity_direction = 2` through the production `-ic` entry path; the
  axis conventions, proof obligations, defects fixed, and known limitations
  are documented in `dev_notes/axial_agnosticity.md`. A permanent z-vertical
  golden-master case is deferred to the terrain-following work. (axial_phase6_straka_walls_smoke)
- Native 2D terrain-following runs: the 2D metric divergence (contravariant/
  J-weighted fluxes), the 2x2 elliptic tensor `J A^T H^-1 A` folded into the
  lap2D cross-term coefficient slots, and a geometry-aware 2D preconditioner
  diagonal. Proven by a 2D div∘correction composition oracle (flat/terrain,
  with/without Coriolis), forced-flat bit-identity, and the native-2D
  resting-atmosphere / conservation gates. (terrain_2d_native)
- Terrain phase 6: advection uses contravariant vertical / J-weighted horizontal
  mass fluxes (HLL upwinding follows automatically), the cell update divides by
  J, recovery's Courant velocity divides J back out, and the CFL gains the
  metric vertical signal speed. Orography is evaluated on periodic-wrapped
  coordinates — a non-periodic hill on a periodic axis otherwise makes the
  elliptic system inconsistent at the duplicated nodes (found via dense
  null-space analysis; bicgstab diverged). smoke_agnesi now runs the 400 m
  witch-of-Agnesi hill end-to-end: J-weighted mass/P conserved to 1e-10,
  max |w| ≈ 0.44 m/s vs the ~0.5 m/s linear estimate. (terrain_phase6_advection_cfl)
- Terrain phase 5: gravity-axis ghost cells respect the terrain — the
  CONTRAVARIANT momentum (mom_v − G·mom_h) is reflected oddly and the
  Cartesian vertical momentum rebuilt with the ghost cell's slope terms;
  ghost hydrostatics use the physical image height and local dz = J·dη.
  Advection sweeps flip the metric alongside the solution arrays. Gates:
  uniform flow over a forced-flat metric matches the plain path to 1e-12. (terrain_phase5_bottom_bc)
- Terrain phase 3: the 3D elliptic operator folds the metric into its tensor —
  `C_ij = wplus ⊙ (J Aᵀ H⁻¹ A)` in role space (lap3D cross terms self-activate)
  and the Helmholtz centre term is J-weighted at nodes. New elliptic oracle
  proves the operator equals the discrete divergence∘momentum-correction
  composition, flat and over an Agnesi hill, to ~1e-12. (terrain_phase3_elliptic_tensor)
- Terrain phase 4: hydrostates at physical height. `States` gains a field mode
  (full per-column fields when terrain is active); `analytical_state` evaluates
  its closed form at z(ξ,η) with local dz = J·dη, `integrated_state` gains a
  fine-grid quadrature branch. Resting-atmosphere oracle: balanced atmosphere
  over a 400 m Agnesi hill stays at rest to ~1e-10 m/s (the solve floor),
  identical to flat. `column`/`initial_pressure` assert no terrain. (terrain_phase4_hydrostates_fields)
- Terrain phase 1: the nodal divergence computes `J∇·F` when terrain is active —
  J-weighted horizontal fluxes and the contravariant vertical flux
  `θ(mom_v − G1·mom_h1 − G2·mom_h2)` in role space, with unchanged differencing
  stencils. Ghost-slab wall zeroing covers the contravariant fluxes automatically.
  Forced-flat oracle (≤1e-13 vs plain path) and role-wiring checks added. (terrain_phase1_metric_divergence)
- Terrain phase 2: pressure gradients map to physical space via the terrain
  gradient matrix A (slope correction of the horizontal rows, 1/J on the
  vertical) in both the explicit forward step and the implicit pressure
  correction; the explicit π update divides the J-weighted divergence by the
  node Jacobian. Bypassed entirely without terrain. (terrain_phase2_metric_gradients)
- Terrain phase 7: Rayleigh sponge generalised — profile builders read the
  configured vertical axis via `axes.coords_along` (η-based taper, correct in
  terrain-following coordinates), damping broadcasts axis-aware for 3D fields,
  the background Y uses field-mode hydrostates under terrain, and the third
  velocity component is damped in 3D. Lifts the documented vertical=1
  limitation; the 2D path (lamb golden masters) is bit-identical. (terrain_phase7_sponge_generalisation)
- General curvilinear advection and consumers (tfc pt 4): the advective mass
  flux is now uniformly rhoY (N_i . m) / rho on every sweep axis (numpy + JAX
  twins share one assembly); CFL uses the contravariant speeds
  |N_a . m|/(rho J) and signal bounds c |N_a|/J on all axes; field-mode
  hydrostates and ghost-cell hydrostatics use the true vertical thickness
  z_eta = J/(N_v)_v instead of J (identical for vertical-line maps). New
  end-to-end Tier-2 gates: resting atmosphere + mountain-wave smoke on an
  x-stretched general map (test_scripts/test_terrain_stretched_smoke.py).
  Legacy G1/G2/z stay stored (consistent with N by construction and asserted
  in tests) — deriving them per sweep flip would cost hot-path allocations. (tfc_phase4_general_advection)
- General curvilinear divergence fluxes (tfc pt 2): the terrain rhs now
  assembles `F_a = N_a . (theta m)` (vertical-first contraction, bit-exact
  reduction to the legacy J-weighted/contravariant fluxes) in both the numpy
  and JAX backends; differencing stencils untouched. New freestream-preservation
  tripwire `test_scripts/test_freestream.py` (uniform flow on a wavy 3D Tier-2
  map: defect bounded and second-order convergent; bit-exact zero on flat). (tfc_phase2_general_divergence)
- General curvilinear gradient map and elliptic fold (tfc pt 3):
  `apply_gradient_map` now applies A_{ka} = (N_a)_k / J and
  `elliptic_tensor(_2d)` folds M = (1/J) N H^-1 N^T (new
  `elliptic_diag_geometric` for the 2D preconditioner diagonal, H^-1 still
  excluded). Bit-exact h == 0 reduction, one-ulp vertical-line reduction;
  laplacian kernels untouched; the JAX device path shares the same functions.
  Oracle suites extended with genuinely stretched Tier-2 maps (3D + native
  2D + Coriolis) and an SPD/symmetry gate. (tfc_phase3_general_elliptic)


Fixed
^^^^^

- Fixed the psinc->comp blending conversion being silently discarded:
  `time_update.do` assigned the pre-conversion `sol`/`npf` aliases returned by
  `prepare_blending` back onto the model state, undoing the conversion's rebind
  (the pre-refactor code threaded the CONVERTED Sol/mpv through). Also pinned
  the throwaway look-ahead steps in `do_psinc_to_comp_conv` /
  `do_lake_to_swe_conv` to the limit-regime clock (`window_step = 0`, matching
  the paper-era `[0, step-1]` call), so under continuous blending the extracted
  half-time pressure is the projected one, not a compressible unprojected one.
  The warm-bubble case's loose 1e-0 momenta tolerances — which had hidden the
  discarded conversion — are tightened to the 1e-5 defaults, and its golden
  master was regenerated deliberately. Visible effect: blended-DA ensembles now
  recover the full balanced pressure field after each assimilation (residual
  acoustic imbalance in p2_nodes is gone). (blending_conversion_discarded)
- Fixed the psinc-to-comp blending conversion being discarded entirely: the
  orchestration bound `sol, npf` aliases BEFORE the conversions and
  `time_update.do` assigned them back over `mem.sol`/`mem.npf`, clobbering the
  reverted-and-converted state that `do_psinc_to_comp_conv` installs. The run
  continued from the throwaway step's half-advanced psinc cells and unblended
  end-of-step pressure instead — the pre-refactor code (`7f0b676~1` schemes.py,
  paper-era `archive/localdab` data.py) threads the CONVERTED `Sol, mpv`
  through. Also runs the throwaway pressure-extraction step on the reference's
  clock (`[0, step-1]`): with the live `window_step == no_of_pi_initial`,
  continuous blending flipped the throwaway to compressible and extracted an
  unprojected pressure — this hit every blended DA analysis window. After the
  fix the blended warm bubble tracks the pure-psinc pressure to 6.6e-5 at the
  blend step (pre-fix 6.6e-4, unblended 3.1e-2). Warm-bubble golden master
  regenerated deliberately (tolerances untouched); fast set, blending-swe and
  blending-hydrostatic verified bit-identical at tol=0. (blending_psinc_to_comp_state_discarded)
- Fixed the comp-psinc blending conversion advancing the real clock: the
  throwaway pressure-extraction step inside `do_psinc_to_comp_conv` reverted
  `mem.sol`/`mem.npf` but not `mem.time`, so every blended step silently
  consumed one dt of integration without advancing the state (pre-ModelState
  code passed t/step by value; `swe_lake.py` already applies the freeze/restore).
  In blended-DA ensembles this made members lag the truth by one dt per
  assimilation window. The warm-bubble golden master was regenerated
  deliberately (the only case exercising this path — fast set, blending-swe and
  blending-hydrostatic remain bit-identical); the inline CompareSol gate passed
  both before and after. (blending_throwaway_clock)
- Clarified the ``CompareSol`` regression-failure message, which read "Relative L2
  error … exceeds tolerance" while the gate is actually per-field max-abs. The gate
  itself (max-abs < tolerance, same 7 fields) is unchanged. (comparesol_message_clarify)
- Fixed the hydrostatic background of the two constant-N stratified regression
  cases (`test_internal_long_wave`, `test_blending_hydrostatic`): their `sol_init`
  now uses the stratification-consistent `hydrostatics.integrated_state` instead of
  the ISOTHERMAL `hydrostatics.analytical_state`. The isothermal background gave
  these constant-N cases an effective N-squared 3.19x too large (gravity-wave
  frequency ~1.79x too fast), a defect masked only because their regression targets
  had been self-generated with the same wrong background. The dominant IGW frequency
  now matches the thesis-era reference to ~1.3%. Targets for both cases regenerated
  (pre-approved bug fix). The five other `alpha_w = 1` gates
  (`test_blending_warm_bubble`, `test_travelling_vortex`, `test_igw_baldauf_brdar`,
  `test_lamb_wave`, `test_unstable_lamb`) do not use the changed function and were
  re-verified bit-identical; `analytical_state` is left untouched (correct for the
  genuinely isothermal cases). See `dev_notes/hydrostatic_blending.md` (ROOT CAUSE). (constant_n_hydrostatic_background)
- Golden-master tolerances recalibrated for cross-platform CI: the first GitHub
  runner pass deviated from locally generated targets by 2.3e-6 (igw rhou) to
  7.0e-5 (Agnesi rhou) — different CPU/BLAS/numba reorder the bicgstab
  reductions, ~100x the same-machine scatter the old gates were tuned to. igw
  returns to the 1e-5 default; the two terrain cases (long elliptic iteration
  chains) gate at 5e-4. Physics remains guarded by the analytic oracles. (cross_platform_tolerances)
- Rebuilt `data_assimilation/prepare.py` on the ModelState/EnsembleState API:
  ensemble members are constructed from fresh `CellSolField`/`NodePressureField`
  containers with one seeded `sol_init` call each (matching the paper-era
  reference, where `sol_init` never pre-ran for N>1 — re-running it on the
  initialised member 0 double-applies the `+=` initialisations and blows up the
  forecast), each member gets its own `FlowSolverCache`, and ghost cells are set
  at construction. Also: observation loading is skipped when `da_times` is empty
  (pure ensemble forecasts need no obs file), `obs_path` now defaults to `None`
  with an actionable error pointing at the dap-rewrite route, and the broken
  `es.flux`/`ensembble_state` references are gone. Verified with an N=2
  travelling-vortex ensemble forecast smoke run. (da_prepare_modelstate)
- Fixed the production rising-bubble input (`-ic rb`): its `sol_init` still
  called the renamed `hydrostatics.state` (now `integrated_state`), so the case
  could not run at all. The MWR-2022 OSSE driver now uses `rb` for the bubble
  experiments (native 160x80, t_ref = 1000 s, seeded delth machinery) instead of
  scaling up `test_blending_warm_bubble`, which is a 31-step blending smoke and
  is CFL-unstable at the paper grid and times. The driver's bubble aux carries
  `CFLfixed` but deliberately not `imbal` — the latter now triggers the initial
  *hydrostatic* conversion from the hydrostatic-blending work, which NaNs this
  nonhydrostatic case; the paper's initial blending is the pseudo-incompressible
  one driven by `initial_blending` alone. (da_rb_input_fix)
- Fixed three run-blockers surfaced by the first N=2 assimilation cycles since
  the DA repair: `dask.diagnostics` is imported explicitly (dask no longer
  auto-imports it), the LETKF's scipy bindings are pinned to the intended
  objects (`scipy.sparse.linalg.spsolve`, dense `scipy.linalg.eigh`,
  `scipy.sparse.eye/diags` — the reference's `import scipy.sparse as sp` only
  resolved under pre-1.8 scipy), and `prepare_rloc` is now built for every
  `da_type` at N>1 because `obs_noiser` needs its cell/node attribute partition
  (the reference had the same latent NameError for ETPF with observation noise).
  LETKF-rloc and ETPF both complete a 2-cycle N=2 travelling-vortex OSSE. (da_runtime_api_fixes)
- Fixed a latent ``TypeError`` in ``EnsembleState.set_members``: the assertion
  ``len(self.set_members == members)`` compared a bound method to a list and then
  called ``len`` on the resulting bool. It now checks ``len(members) == len(self.members)``
  (the ensemble size is preserved when forecast members are replaced by analysis members). (ensemblestate_set_members_typeerror)
- Fixed the incompressible + field-mode gravity ghost fill in
  `cell_boundary._calculate_ghost_values`: it indexed the field-mode
  `HydroState.rhoY0` (and, in the `ATMOSPHERIC_EXTENSION` branch, `p20`) with
  only the vertical component `nimage[y_axs]`. That is correct for the 1D
  PROFILE-mode hydrostates (vertical profiles broadcast on demand) but wrong
  for the full grid-shaped FIELD-mode hydrostates a terrain / sphere run
  carries — a scalar radial index then slices the wrong (longitude) axis and
  mis-shapes the result, raising on the `rho = rhoY * S` broadcast (or filling
  ghosts with garbage where the axes happen to match). The branch runs only
  during `do_initial_projection` (which freezes the regime to incompressible),
  so no prior case hit it: the SWE sphere cases that project have `grav = 0`,
  so the gravity ghost fill never runs. Now indexes the whole ghost slice
  tuple when the hydrostate is in field mode (via a `_hydro_at` helper), like
  the sibling `metric.height[idx]` / `sol.rhoY[idx]` branches. Regression:
  `test_scripts/test_field_mode_gravity_ghost.py`. This unblocks the
  TC2-style initial projection on the deep compressible shell (Hughes &
  Jablonowski pt 2). (field_mode_gravity_ghost)
- Fixed the `rb` initial-condition registry entry (missing `pybella.` package prefix) and removed 18 dead legacy keys from `IC_MODULES` that pointed at modules deleted in the package restructure; the legacy ICs remain recoverable from the `archive/full_coriolis` tag. (ic_registry_cleanup)
- Fixed a latent axial-agnosticity violation in the Lamb-wave regression cases:
  `test_lamb_wave` / `test_unstable_lamb` hard-coded `ud.bdry_type[1]` when applying
  the Rayleigh sponge boundary instead of `ud.bdry_type[axes.vertical_axis(ud)]`
  (inert today as both run vertical=1, wrong for any other gravity direction). The
  four Rayleigh-switch sites now route through a single
  `tests/case_setup.apply_rayleigh_bdry` helper (axis from `axes.vertical_axis`),
  with `with_tau` folding in the terrain cases' `get_tau_y`. The iny-resizing
  `rayleigh_bc_function` stays vertical=1 (documented limitation). Bit-identical at
  tol 0 on all four cases. (rayleigh_bdry_helper_axis)
- CI: drop the redundant `test_jax_*` pytest invocations from the numpy
  integration job — without jax they collect zero tests, which pytest >= 9
  fails with exit code 5; the jax-equivalence job runs them. (ci_jax_pytest_exit5)
- Fixed the 2D out-of-plane Coriolis defect: the `rhow` momentum row in
  `explicit_euler.do_forward_step` was guarded by `if ndim == 3`, so 2D runs
  applied only the implicit half of the out-of-plane Coriolis rotation. Found
  by the Baldauf-Brdar analytic oracle (out-of-plane velocity error pinned at
  ~0.44 rel-L2 independent of dt, sim/ref amplitude ratio ~0.6, with O(f t)
  feedback into u); after the fix the error drops to 0.06 (dt-convergent to
  0.03) and the amplitude ratio to 1.03. Affects only 2D runs with Coriolis
  components in the x/y slots: the `igw_baldauf_brdar` and
  `internal_long_wave` golden-master targets were deliberately regenerated;
  all other cases are bit-identical (Lamb cases use only `strength[2]`, which
  does not enter the w-row). Oracle gates tightened accordingly
  (`test_igw_analytic.GATES`: vo 0.60 -> 0.10). (fix_2d_coriolis_w_row)
- Fixed the full-3D (`inz > 1`) implicit/elliptic solver path, broken since the
  package restructure: re-wired `lap3D` to the interior-sized `npf` arrays
  (missing imports, coefficient slicing, 3D preconditioner via
  `preconditioner.prepare_diag`), fixed the flat-vector memory layout
  (C-order `[x, y, z]`; the old reshape silently transposed x and z on
  non-cubic grids), fixed the 3D `rhs` shape mismatch and a sign error on the
  x-component of the 3D nodal divergence (legacy, dating to Oct 2021), and made
  the pressure diagnostic kernel dimension-agnostic. The 3D path is validated
  against the 2D solver on a y-uniform quasi-2D problem to ~1e-10
  (`test_scripts/test_3d_elliptic_oracle.py`), and the
  `test_travelling_vortex_3d_coriolis` regression case now runs with a
  committed golden-master target. (fix_3d_elliptic_path)
- Fixed a latent view bug in the `lap3D` periodic ghost reconstruction: the
  tmp-swap of the duplicated periodic rows used a numba *view* (`tmp = p[1]`),
  making the closing write-back a no-op — duplicates were mirrored instead of
  exchanged and the operator column at row 1 was dead. Proven outcome-inert
  (end-to-end solutions bit-identical: the solver only visits periodically
  consistent vectors, where the two conventions coincide exactly), so no
  regression targets move. The numpy kernel (`tmp = p[1].copy()`) and the JAX
  twin's permutation were fixed in lockstep; equivalence suite and all 3D
  oracles pass. (lap3d_periodic_view_swap)
- `lap2D_manual` now treats RAYLEIGH boundaries as walls, matching `lap3D`,
  `lap2D_numba` and the divergence slab-zeroing. Previously a sponged top in a
  native-2D run got a periodic-in-y elliptic stencil. No-op for the existing 2D
  golden masters (their vertical handling goes through the atmospheric-extension
  branch — verified bit-identical, and the stable lamb wave still propagates at
  0.99 Cs with stable amplitude at twice the regression horizon). (lap2d_rayleigh_wall)


Improved Documentation
^^^^^^^^^^^^^^^^^^^^^^

- Repo-root hygiene: stale SHA-migration scratch removed; the Phase-C GPU report and
  TFC generalization plan moved under `dev_notes/`; doc references updated. (root_scratch_hygiene)
- Added ``STYLE_GUIDE.md`` — a committed, human- and AI-readable coding guide
  covering the hard invariants, the reproducibility-gate "iron rule" and gate
  tiers (``test_scripts/repro_gate.py``), structure/hygiene/test conventions, and
  an agent validation checklist (which gate proves which kind of change). (style_guide)


0.60.0 (2024-12-23)
-------------------

Changed
^^^^^^^

- Restructured outermost looping in main() (23a7817db225fe2188bcf3045cfdc863435b6888)


0.50.6 (2024-12-10)
-------------------

Changed
^^^^^^^

- Cleaned __main__ with restructuring of code (61e48c2c60e32aa8ac7699979a9836c8d17ded64)


0.50.5 (2024-03-24)
-------------------

Infrastructure
^^^^^^^^^^^^^^

- Standardised import statements throughout codebase (#25)


Improved Documentation
^^^^^^^^^^^^^^^^^^^^^^

- Reviewed and updated Readme (#10, #17)
- Reincluded Sphinx docs (#23)


0.50.1 (2024-03-21)
-------------------

Improved Documentation
^^^^^^^^^^^^^^^^^^^^^^

- moved logging initialisation to io.py; initialiser auto create directory; log filename now includes run case and datetime stamp (#19)
- added changelog function (#22)
