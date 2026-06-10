import pytest
import subprocess


@pytest.mark.parametrize(
    "ic",
    [
        "test_travelling_vortex",
        "test_internal_long_wave",
        "test_igw_baldauf_brdar",
        "test_lamb_wave",
        "test_unstable_lamb",
        "test_swe_vortex",
    ],
)
def test_single_run(ic):
    result = subprocess.run(
        ["pybella", "-ic", ic, "-N", "1"], capture_output=True, text=True
    )

    assert result.returncode == 0, (
        f"Command failed with return code {result.returncode}\n"
        f"STDERR:\n{result.stderr.strip()}\n"
        f"STDOUT:\n{result.stdout.strip()}"
    )
