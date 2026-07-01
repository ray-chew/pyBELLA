Reinstated nonhydrostatic<->hydrostatic dynamics-regime blending (thesis ch. 4.2)
as a schedule-only mechanism, with no explicit conversion routine (matching the
thesis-era RKLM_Python). The eos schedule (`physics/eos.py`) flips the blending
parameter `alpha_w` (`is_nonhydrostatic`/`nonhydrostasy`) per step -- hydrostatic
for the first `no_of_hy_initial` steps, then nonhydrostatic -- and `alpha_w` is
wired into the two places that make the vertical-momentum equation degenerate to
hydrostatic balance when it is zero: the elliptic vertical self-coupling
`h22 = 1/nu_nh` (`numerics/coriolis.py` + JAX twin; drops the spurious `nonhydro`
factor, giving the thesis coefficient of eq. 4.40/4.42) and the explicit
vertical-momentum inertia (`numerics/explicit_euler.py`). The same
predictor->corrector step runs every timestep; `alpha_w` alone selects the regime.
Added the `test_blending_hydrostatic` regression case (imbalanced Skamarock-Klemp
hydrostatic inertia-gravity wave) and wired it into CI. Reproduces the thesis
sec. 6.2.3 w-probe suppression (centre-probe vertical-velocity rel-err O(0.01-0.1),
matching the reported 0.024 within the harness normalisation). Bit-identical for
`alpha_w = 1`. See `dev_notes/hydrostatic_blending.md`.
