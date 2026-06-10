Fixed the 2D out-of-plane Coriolis defect: the `rhow` momentum row in
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
(`test_igw_analytic.GATES`: vo 0.60 -> 0.10).
