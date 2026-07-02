SWE <-> lake blending reinstated on the ModelState API (2D x-y layout):
`do_swe_to_lake_conv` / `do_lake_to_swe_conv` ported, the orchestration's
broken SWE branches repaired (undefined `flux`, unbound conversion flag,
duplicate-`ud` latent TypeError on the continuous-blending call), and a new
golden-master regression case `test_blending_swe` added (balanced SWE vortex
with initial swe->lake->swe blend).
