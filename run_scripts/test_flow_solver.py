import pytest
import subprocess

@pytest.mark.parametrize("ic", 
                        ["test_travelling_vortex",
                        "test_internal_long_wave",
                        "test_lamb_wave"])
def test_single_run(ic):
    result = subprocess.run(
    ["pybella", "-ic", ic, "-N", "1"],
    capture_output=True,
    text=True
    )

    assert result.returncode == 0, result.stderr.splitlines()[-3:]