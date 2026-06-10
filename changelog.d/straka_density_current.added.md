Added the Straka density current (Straka et al. 1993) as regression case
`test_straka` — the suite's first nonlinear, advection-dominated gravity case
(256x32 at 200 m, dt = 4 s to t = 900 s; front position, peak winds and
symmetry verified against the published benchmark). Includes a new explicit
constant-coefficient diffusion module (`flow_solver/numerics/diffusion.py`,
enabled per-case via `ud.diffusion` / `ud.diffusion_coeff`, off by default)
implementing the benchmark's fixed K = 75 m^2/s on velocity and potential
temperature. The case runs periodic in x: the x-WALL elliptic path is
currently broken/untested (wall momentum zeroing in the nodal divergence
handles the vertical axis only) — known limitation, to be fixed with the
axial-agnosticity refactor.
