Added ``test_scripts/repro_gate.py``, a reproducibility-gate harness for
behaviour-preserving refactors: ``capture`` snapshots each case's output H5 and
per-field max-abs to a scratch baseline, ``check`` reruns and requires the inline
CompareSol gate to stay green *and* the output to be bit-identical (tol 0) to the
baseline. Bit-identity is the strong oracle proving a refactor changed no numerics.
