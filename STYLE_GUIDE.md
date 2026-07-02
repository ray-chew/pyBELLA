# pyBELLA — Coding Style & Contribution Guide

> Human- and AI-readable. If you are an automated agent: this file is normative.
> Read the **Hard invariants** and **The iron rule** before touching numerics.
> Cross-references: project setup in [CLAUDE.md](CLAUDE.md); design notes in
> [dev_notes/](dev_notes/); a navigable code graph in `graphify-out/`.

pyBELLA is a compressible atmospheric-flow solver (blended dynamics + ensemble
DA) in pure numpy + numba, with a parallel JAX backend, validated bit-for-bit
against numpy golden masters. These conventions exist to keep that validation
true and the hot path fast.

---

## Hard invariants (a change that violates one is wrong, not just unstylish)

1. **numba hot path.** ~28 `@njit`/`@stencil` kernels (`utils/operators/*`,
   `numerics/explicit_advection/*`, `numerics/coriolis.py`,
   `physics/hydrostatics.py`). No Python objects/classes/dicts/closures inside
   kernels; no dynamic typing. Extraction is fine **if** the extracted helper
   stays njit-compatible (numba inlines small njit functions).
2. **In-place mutation is the norm.** Solver state is updated via `arr[...] =`,
   `+=` on `ModelState` fields. Do not "functionalise" the numpy path; that is
   the JAX backend's job, and the asymmetry is structural.
3. **Elliptic solve = scipy `bicgstab` + diagonal preconditioner**
   (`numerics/implicit_euler.py`, `utils/operators/laplacian/preconditioner.py`).
   The hardest piece to change; touch with a parity test in hand.
4. **numpy↔JAX bit-parity.** `backends/jax_ops/*` mirrors the numpy core and is
   validated bit-for-bit (H100). Duplication between the two is often
   *intentional mirroring* — a shared abstraction is only acceptable if it
   cannot perturb parity. Never fix a bug on one side only.
5. **Axial agnosticity.** `utils/axes.py` `role_perm(v)` is the single source of
   truth for the `(h1, v, h2)` role-space ↔ array-axis mapping (cyclic perms
   only — Omega is a pseudovector). Do not hardcode axis indices.
6. **Terrain reduction contract.** `elem.metric is None` ⇒ the uniform-Cartesian
   path is **bit-identical**. Metric terms flow through the J-weighted
   divergence, A-mapped gradients, and the elliptic tensor only when present.
7. **TFC curvilinear contract.** The Klein curvilinear generalization keeps
   legacy `G1/G2/z` arrays stored intentionally; the general path must reduce
   bit-exactly to the legacy formulas for a vertical-line metric. See
   `dev_notes/tfc_generalization_plan.md`.

---

## The iron rule: the reproducibility gates are the oracle

A change that is meant to change **no numerics** must leave every solver output
**bit-identical**. This is proven, never asserted, and never worked around.

- **Never** to make a refactor "pass": regenerate a target (`diag_updt_targets`
  stays `False`), loosen a `tolerance`, or change which fields `CompareSol`
  compares. Those move the pass criteria and hide regressions.
- The inline gate (`tests/diagnostics.py::CompareSol.test_do`) raises when a
  per-field max-abs error meets/exceeds its tolerance, making `pybella` exit
  nonzero. **Exit 0 == gate green.** The 7 gated fields are
  `rho, rhou, rhov, rhow, rhoY, rhoX, p2_nodes`; the real reference is the
  stripped HDF5, not `test_targets.yml`.

### Gate tiers and the anti-tamper ritual

Use `test_scripts/repro_gate.py` (the protocol in code):

```bash
python test_scripts/repro_gate.py capture --set fast   # before a refactor
# ... make the change ...
python test_scripts/repro_gate.py check   --set fast   # exit 0 == inert
```

- **fast** — `travelling_vortex`, `straka`, `travelling_vortex_3d_coriolis`
  (2D-elliptic / diffusion+walls / 3D+Coriolis). Cheap; run after every commit.
- **terrain** — `agnesi_hydrostatic`, `schaer_ridge`. Slow (per-process numba
  JIT of the 3D terrain kernels); run at phase boundaries.
- **full** — every golden-master case; pre-merge.
- **parity** (backend changes only) — rerun under `PYBELLA_BACKEND=jax` and
  `jax-device` (`JAX_ENABLE_X64=1`) + the `test_jax_*` equivalence/fullrun
  scripts.

The strong oracle is **bit-identity at tol 0** (`compare_h5_runs.py`): a
byte-identical output H5 implies the gate verdict and every max-abs are
unchanged. For `jax-device` (documented Krylov floor) parity to numpy is *not*
bit-exact; the gate is CompareSol < 1e-5 vs the same targets, plus
device↔device bit-identity (baseline vs candidate) to prove a split mechanical.

When the change has no numerical output to diff (e.g. a value-preserving data
refactor), add a second oracle that *does* prove value-identity — e.g. the
`vars(UserData())` attribute-table diff used for the `case_setup` refactor.

---

## Structure conventions

- **No god-objects / god-modules.** When a module/class accretes unrelated
  responsibilities, split it into a package of cohesive submodules and re-export
  the public names from `__init__` so call sites are unchanged. Worked example:
  `utils/io/` (`writer` / `restart` / `debug` / `cli`), split from the former
  800-line `io.py` with a re-exporting `__init__`.
- **Mechanical splits stay mechanical.** Move code verbatim (extract exact line
  ranges; don't retype), fix only relative-import depth, and prove inertness
  with the bit-diff gate. A split that "improves" code while moving it is two
  changes — do them separately.
- **Backend dispatch goes through the helpers.** Detect the backend with
  `backends.is_jax_backend(ud)` / `is_device_backend(ud)` — never inline
  `getattr(ud, "backend", "numpy") in (...)`. The helpers are the single place
  the backend-name set is defined.
- **The `jax_ops` mirror.** Every numpy operator with a JAX twin keeps the same
  public signature; structural reorganisation of the JAX side must preserve the
  call sequence (and the device compile-cache) so parity holds.

---

## Local hygiene

- **Named constants over magic numbers.** A literal with physical or historical
  meaning gets a named module constant with a one-line comment (e.g.
  `_CFLFIXED_DT_SECONDS = 21.69`).
- **Delete dead code; don't comment it out.** Confirm zero live references
  first (`grep`), then remove. History lives in git.
- **Assertions must be able to fire.** `assert len(self.foo == bar)` (a bound
  method compared to a value) is a bug, not a check.
- **Format with `black`** before committing. No repo-wide reformatting in an
  otherwise-behavioural change — it bloats the diff and defeats tol-0 review.

---

## Test conventions

- **A regression case** is `tests/test_<case>.py` with a `UserData` class (no
  args; consumed as `vars(UserData())`) and `sol_init(Sol, npf, elem, node, th,
  ud, seed=None)`, registered in `interfaces/ic_config.py::IC_MODULES`. It must
  appear in CI (`test_scripts/test_flow_solver.py` or `test_blending.py`).
- **Case physics is per-case; only mechanical idioms are shared.** Use
  `tests/case_setup.py`: `build_bdry(...)` (the per-instance, object-dtype
  boundary triple — never a shared class attribute) and `make_diag_state(...)`
  (centralises the `Nx=inx-1`/`Ny=iny-1`/`steps=[stepmax-1]` offsets, forwards
  case-specific kwargs). Do **not** try to share `stratification_function` /
  `rhoe_function` — they encode case-specific physics.
- **(Re)generating a target** is a deliberate, separate act (never part of a
  refactor): set `diag_updt_targets = True`, run the case, commit the stripped
  H5 + PNGs, set it back to `False`, re-run to confirm the gate asserts.
- **`time_increment=True`** cases compare the *difference* between consecutive
  steps for pressure fields, not absolute fields.

---

## Tooling & workflow

- **Branch** off `develop`; feature branches `ext/*`, `optimise/*`, etc.
- **One towncrier fragment per change** in `changelog.d/<slug>.<type>.md`
  (`added`/`changed`/`fixed`/`removed`/`deprecated`/`docs`/`infrastructure`).
  `changelog.d/` is gitignored, so add fragments with `git add -f`.
- **Commit style:** short title + 3–6 bullets, no attribution trailers.

### AI-agent validation checklist (which gate for which change)

| Change kind | Required proof |
|---|---|
| value-preserving data refactor (e.g. UserData) | attribute-table identical **+** FAST bit-diff |
| numpy structural refactor (split/move, no math) | FAST bit-diff tol 0 (+ FULL at phase boundary) |
| anything touching `backends/` or dispatch | numpy bit-diff **+** jax & jax-device green (+ device↔device tol 0 for a device split) |
| numerics change (intended behaviour change) | new target generated deliberately; reviewed; FULL gate |
| dead-code removal | re-grep zero live refs **+** FAST bit-diff tol 0 |
