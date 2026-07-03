"""Native ensemble data assimilation: LETKF (batch / r-localised) and ETPF.

FROZEN at the Chew-Benacchio-Klein 2022 (MWR) reproduction (2026-07):
2D x-y (vertical = 1), numpy backend only, bug fixes only. The maintained DA
engine is the NEDAS interface (dev_notes/nedas_interface.md), which replaces
this layer; the OSSE parity targets are pinned in
dev_notes/da_reinstatement.md and regenerable via run_scripts/osse_mwr2022.py.
Regression guard: test_scripts/test_da_smoke.py (CI).
"""
