Reinstated nonhydrostatic<->hydrostatic dynamics-regime blending (thesis ch. 4.2)
as a schedule-only mechanism, with no explicit conversion routine (matching the
thesis-era RKLM_Python). The eos schedule (`physics/eos.py`) flips the blending
parameter `alpha_w` (`is_nonhydrostatic`/`nonhydrostasy`) per step -- hydrostatic
for the first `no_of_hy_initial` steps, then nonhydrostatic -- and `alpha_w` is
wired into the two places that make the vertical-momentum equation degenerate to
hydrostatic balance when it is zero: the elliptic vertical self-coupling
`h22 = 1/nu_nh` (`numerics/coriolis.py` + JAX twin; drops the spurious `nonhydro`
factor, giving the thesis coefficient of eq. 4.40/4.42) and the explicit
vertical-momentum inertia (`numerics/explicit_euler.py`). Nonhydrostatic steps
(`alpha_w = 1`) run the normal predictor->corrector; the hydrostatic step
(`alpha_w = 0`) runs the thesis one-step recipe (sec. 4.2.3-4.2.4): two first-order
predictor updates for an input-imbalance-independent balanced vertical momentum
plus a second-order `pi` update, then the second-order corrector.
Added the `test_blending_hydrostatic` regression case (imbalanced Skamarock-Klemp
hydrostatic inertia-gravity wave) and wired it into CI. Reproduces the thesis
sec. 6.2.3 w-probe suppression (centre-probe vertical-velocity rel-err ~0.022,
matching the reported 0.024 within the harness normalisation). Bit-identical for
`alpha_w = 1`. See `dev_notes/hydrostatic_blending.md`.
