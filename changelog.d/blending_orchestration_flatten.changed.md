Flattened the blending orchestration: named predicate helpers
(`_window_start_conversion_due`, `_full_blend_due`, `_initial_blend_phase`)
replace the nested conditionals, the field-order-dependent `ModelState`
iterator unpack is gone, and the dead `debug` threading is dropped.
Bit-identical on all three blending masters and the fast set.
