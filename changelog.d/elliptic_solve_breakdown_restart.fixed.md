Elliptic solve: silent bicgstab breakdown fixed. scipy's ``info`` return was
discarded, so a breakdown (``info < 0``) or a maxiter exit handed the
unconverged iterate to the momentum correction as if it were the solution. On
the sphere the pure-Neumann initial projection is ill-conditioned enough (a
cluster of near-null modes from the 1/cos(phi) metric; condition 1e5-1e6 on
the coarse grids vs 1e3 Cartesian) that BiCGSTAB broke down after a few
hundred iterations on every sphere SWE case, at relative residuals of 1e-4 to
5e-2 that depended on the BLAS thread count — all four sphere golden masters
embedded a failed projection, and CI (single-threaded runner) rightly
disagreed with them by up to 1.4e-3. The solve now restarts from the current
iterate with a fresh shadow residual on breakdown (one restart converges every
case; thread-1 and thread-24 answers then agree to ~1e-7), the achieved
residual is checked after every solve (raise on breakdown, warn on maxiter,
within a 10x slack of the stopping criterion), ``ud.rtol`` exposes the
relative tolerance that actually governs convergence (scipy default 1e-5,
unchanged), and the hybrid JAX solve mirrors the restart and residual check.
The four sphere SWE golden masters are regenerated from converged projections.
