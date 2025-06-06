import pytest
import subprocess

@pytest.mark.parametrize(
    "ic",
    [
        "test_blending_warm_bubble",
    ],
)
def test_single_run(ic):
    result = subprocess.run(
        ["pybella", "-ic", ic, "-N", "1"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr.splitlines()[-3:]
