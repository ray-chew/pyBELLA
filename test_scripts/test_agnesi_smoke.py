"""Plumbing smoke test: quasi-2D mountain-wave configuration shape.

Runs the ``smoke_agnesi`` case (x horizontal, y vertical, z degenerate
periodic) via the real ``pybella -ic`` entry point for a few steps —
pushing a vertical-slice setup through the full-tensor 3D elliptic path
(``lap3D``) that terrain metric terms attach to. Flat in Phase 0 of the
terrain work; the Agnesi hill switches on once the metric-aware operators
land. No golden master: the Agnesi regression case is separate.
"""

import subprocess


def test_agnesi_smoke_production_path():
    result = subprocess.run(
        ["pybella", "-ic", "smoke_agnesi", "-N", "1"], capture_output=True, text=True
    )
    assert result.returncode == 0, (
        f"Command failed with return code {result.returncode}\n"
        f"STDERR:\n{result.stderr.strip()}\n"
        f"STDOUT:\n{result.stdout.strip()}"
    )
