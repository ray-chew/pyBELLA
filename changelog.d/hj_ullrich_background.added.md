Hughes & Jablonowski (2023) mountain-induced baroclinic wave (sphere
follow-up), foundation: the Ullrich et al. (2014/2016) dry balanced base
state as a reusable analytic module (`tests/ullrich_baroclinic.py`) — the
closed-form temperature, pressure, density and gradient-wind zonal jet
T/u/p/rho(phi, z) (Appendix B, dry variant q_v=0). pyBELLA being
height-based, these evaluate directly with no pressure root-finding. The
module also carries the two-ridge H&J topography (Eq. 1, Table 1) with
analytic longitude/latitude gradients (required by `SphericalTerrainMap`).
Gated by `test_hj_background.py`: surface T endpoints (310/240 K), the
~28 m/s midlatitude jet near 45 deg / 10 km, hydrostatic (dp/dz = -rho g)
and gradient-wind self-consistency; the 2000 m ridge peaks, their 10%
half-widths, longitude periodicity and finite-difference-checked gradients.
Next: the field-mode hydrostatic reference + well-balanced sphere IC, the
ridges wired through `SphericalTerrainMap`, and the multi-day run.
