"""Plumbing smoke test: gravity_direction = 2 through the production path.

Runs the ``smoke_zvert`` case (z-vertical, met convention) via the real
``pybella -ic`` entry point for a few steps — covering the IC registry,
``prepare.initialise`` (incl. ``axes.validate``), hydrostatics along z,
axis-2 gravity ghost cells / wall zeroing, and the full time loop, which
the in-process permutation oracle bypasses. No golden master: physics
agnosticity is proven by ``test_permutation_oracle.py``; there is no
permanent z-vertical regression case.
"""

import subprocess


def test_zvert_production_path():
    result = subprocess.run(
        ["pybella", "-ic", "smoke_zvert", "-N", "1"], capture_output=True, text=True
    )
    assert result.returncode == 0, (
        f"Command failed with return code {result.returncode}\n"
        f"STDERR:\n{result.stderr.strip()}\n"
        f"STDOUT:\n{result.stdout.strip()}"
    )
