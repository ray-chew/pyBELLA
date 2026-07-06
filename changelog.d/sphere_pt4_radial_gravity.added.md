Spherical geometry, stage D: radial gravity on the true-radius 3D shell —
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
projection within its solver-floor equivalence class).
