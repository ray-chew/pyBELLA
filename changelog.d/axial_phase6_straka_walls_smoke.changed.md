Axial-agnosticity Phase 6 wrap-up: Straka now runs its faithful free-slip
wall configuration on all boundaries, exercising the repaired x-WALL
elliptic path (target regenerated for the wall config); a z-vertical
plumbing smoke case (`smoke_zvert` + `test_scripts/test_zvert_smoke.py`)
pushes `gravity_direction = 2` through the production `-ic` entry path; the
axis conventions, proof obligations, defects fixed, and known limitations
are documented in `dev_notes/axial_agnosticity.md`. A permanent z-vertical
golden-master case is deferred to the terrain-following work.
