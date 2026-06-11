Fixed a latent view bug in the `lap3D` periodic ghost reconstruction: the
tmp-swap of the duplicated periodic rows used a numba *view* (`tmp = p[1]`),
making the closing write-back a no-op — duplicates were mirrored instead of
exchanged and the operator column at row 1 was dead. Proven outcome-inert
(end-to-end solutions bit-identical: the solver only visits periodically
consistent vectors, where the two conventions coincide exactly), so no
regression targets move. The numpy kernel (`tmp = p[1].copy()`) and the JAX
twin's permutation were fixed in lockstep; equivalence suite and all 3D
oracles pass.
