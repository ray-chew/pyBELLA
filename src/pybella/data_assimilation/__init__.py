"""Native ensemble data assimilation: LETKF (batch / r-localised) and ETPF.

FROZEN at the Chew-Benacchio-Klein 2022 (MWR) reproduction: 2D x-y
(vertical = 1), numpy backend only, bug fixes only, do not extend. The
maintained DA engine is the NEDAS interface (``pybella.interfaces.nedas``),
which replaces this layer; the OSSE parity targets are regenerable via
run_scripts/osse_mwr2022.py.
Regression guard: test_scripts/test_da_smoke.py (CI).
"""
