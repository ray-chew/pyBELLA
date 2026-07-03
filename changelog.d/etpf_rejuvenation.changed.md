Reinstated the ETPF rejuvenation term (`+ delta * randn(...)` after the
transport step), which was commented out in the reference code and caused
severe ensemble collapse (spread ~5e-4 at rejuvenation_factor 0.001) in the
MWR-2022 OSSE reproduction.
