Generalized terrain metric (tfc pt 1): `MetricFields` now carries the Klein
area normals `N` (outer index = array axis, rotating with sweep flips; inner
index = fixed Cartesian component) and physical coordinates `x`, synthesized
bit-exactly from the vertical-line map when not supplied. New `CurvilinearMap`
base + `build_metric_fields_from_map` general-path builder (Tier-2 capable,
J > 0 guarded, effective slopes for legacy consumers); device config mirrors
the normals.
