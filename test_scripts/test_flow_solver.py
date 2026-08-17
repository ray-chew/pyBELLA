import pytest
import subprocess


@pytest.mark.parametrize(
    "ic",
    [
        "test_travelling_vortex",
        "test_travelling_vortex_3d_coriolis",
        "test_internal_long_wave",
        "test_igw_baldauf_brdar",
        "test_lamb_wave",
        "test_unstable_lamb",
        "test_swe_vortex",
        "test_sphere_swe_tc2",
        "test_sphere_swe_tc6",
        # pole-to-pole (BdryType.POLE + polar filter + elliptic ring collapse).
        # The JAX jobs already gate these for jax-vs-numpy EQUIVALENCE; these
        # entries gate the numpy physics against a stored reference, so a numpy
        # pole regression cannot pass silently.
        "test_sphere_swe_tc2_global",
        "test_sphere_swe_tc6_global",
        "test_straka",
        "test_agnesi_hydrostatic",
        "test_schaer_ridge",
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
