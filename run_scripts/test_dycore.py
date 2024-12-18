import pytest
import sys
import subprocess

@pytest.mark.parametrize("ic", 
                        ["test_travelling_vortex",
                        "test_internal_long_wave",
                        "test_lamb_wave"])
# def test_single_run(ic):
    # run = subprocess.Popen(
    #     [sys.executable, "src", "-ic", ic, "-N", "1"],
    #     stdout=subprocess.PIPE,
    #     stderr=subprocess.PIPE,
    # )
    # _, stderr = run.communicate()
    # assert run.returncode == 0, stderr.splitlines()[-3:]

def test_single_run(ic):
    result = subprocess.run(
    ["pybella", "-ic", ic, "-N", "1"],
    capture_output=True,
    text=True
    )

    assert result.returncode == 0, result.stderr.splitlines()[-3:]