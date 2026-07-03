Fixed the comp-psinc blending conversion advancing the real clock: the
throwaway pressure-extraction step inside `do_psinc_to_comp_conv` reverted
`mem.sol`/`mem.npf` but not `mem.time`, so every blended step silently
consumed one dt of integration without advancing the state (pre-ModelState
code passed t/step by value; `swe_lake.py` already applies the freeze/restore).
In blended-DA ensembles this made members lag the truth by one dt per
assimilation window. The warm-bubble golden master was regenerated
deliberately (the only case exercising this path — fast set, blending-swe and
blending-hydrostatic remain bit-identical); the inline CompareSol gate passed
both before and after.
