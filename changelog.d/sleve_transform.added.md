SLEVE vertical transform (Schär et al. 2002; Leuenberger et al. 2010
exponent n): two-scale orography split `ud.orography_smooth` (large-scale)
+ residual against the total `ud.orography`, per-component sinh decay,
analytic Jacobian — the first eta-dependent J through the operators.
Selected via `ud.vertical_transform = SLEVETransform(s1, s2, n)`; the
activation contract and the Gal-Chen path are untouched. Proven by
transform identities, scale-separation property, and SLEVE-parametrized
elliptic-composition + resting-atmosphere oracles (3D and native 2D).
