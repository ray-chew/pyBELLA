import pytest
import sys
import subprocess

@pytest.mark.parametrize("ic", 
                        ["test_travelling_vortex",
                        "test_internal_long_wave",
                        "test_lamb_wave"])
def test_single_run(ic):
    run = subprocess.Popen(
        [sys.executable, "src/__main__.py", "-ic", ic, "-N", "1"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    _, stderr = run.communicate()
    assert run.returncode == 0, stderr.splitlines()[-3:]