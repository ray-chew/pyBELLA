"""Plumbing smoke test: Schär ridge through the production entry point.

Runs ``test_schaer_ridge`` (native 2D, SLEVE transform, RAYLEIGH top) via
``pybella -ic`` — the same invocation the golden-master comparison uses —
and gates on a clean return code. The physics gates live in
``test_schaer_analytic.py``.
"""

import subprocess


def test_schaer_smoke_production_path():
    result = subprocess.run(
        ["pybella", "-ic", "test_schaer_ridge", "-N", "1"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"Command failed with return code {result.returncode}\n"
        f"STDERR:\n{result.stderr.strip()}\n"
        f"STDOUT:\n{result.stdout.strip()}"
    )
