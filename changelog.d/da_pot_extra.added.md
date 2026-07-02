The ETPF's optimal-transport step now uses POT (`ot.emd`) instead of the
undeclared, unmaintained `pyemd` dependency, installable via
`pip install "pybella[da]"` and import-guarded so the deterministic solver and
the LETKF never need it. Both are exact solvers: on an ETPF-shaped random
N=10 problem the transport plans agree to machine epsilon (max abs difference
2.8e-17, identical objective cost), and `emd_with_flow`'s extra-mass penalty
was irrelevant because both ETPF histograms sum to one.
