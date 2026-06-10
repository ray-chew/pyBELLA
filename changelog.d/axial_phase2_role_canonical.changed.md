Axial-agnosticity Phase 2: Coriolis and explicit dynamics are now
role-canonical — `multiply_inverse_terms` binds (wh1, wv, wh2) and the
momentum components in (h1, vertical, h2) role order via the cyclic axis
permutation (njit kernels unchanged); new `compute_inverse_coefficients`
exposes the cached H^-1 fields for the upcoming tensor elliptic operator;
buoyancy and the rhoX stratification coupling act on the configured vertical
momentum; the explicit momentum rows are written in role symbols with
expression trees preserved; the advection pwchi special case keys on the
vertical axis. Bit-identical for all existing cases (verified at tol=0).
