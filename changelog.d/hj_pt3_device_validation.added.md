Hughes & Jablonowski (2023) mountain baroclinic wave, pt 3 (device path):
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
qualitative Rossby wave train vs the paper's dry FV/SE panels + self-convergence.
