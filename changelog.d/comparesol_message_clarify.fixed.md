Clarified the ``CompareSol`` regression-failure message, which read "Relative L2
error … exceeds tolerance" while the gate is actually per-field max-abs. The gate
itself (max-abs < tolerance, same 7 fields) is unchanged.
