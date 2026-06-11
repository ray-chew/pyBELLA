Golden-master tolerances recalibrated for cross-platform CI: the first GitHub
runner pass deviated from locally generated targets by 2.3e-6 (igw rhou) to
7.0e-5 (Agnesi rhou) — different CPU/BLAS/numba reorder the bicgstab
reductions, ~100x the same-machine scatter the old gates were tuned to. igw
returns to the 1e-5 default; the two terrain cases (long elliptic iteration
chains) gate at 5e-4. Physics remains guarded by the analytic oracles.
