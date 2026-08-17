# Initial-condition registry: -ic keys -> module paths providing UserData + sol_init.
#
# ICs from the original RKLM_Python/inputs tree (acoustic waves, SWE variants,
# baroclinic instability, ...) were removed when that tree was dropped; they remain
# recoverable from the git tag `archive/full_coriolis` and are mostly small tweaks
# of the cases below.
IC_MODULES = {
    "rb": "pybella.inputs.rising_bubble",
    "test_travelling_vortex": "pybella.tests.test_travelling_vortex",
    "test_travelling_vortex_3d_coriolis": "pybella.tests.test_travelling_vortex_3d_coriolis",
    "test_internal_long_wave": "pybella.tests.test_internal_long_wave",
    "test_igw_baldauf_brdar": "pybella.tests.test_igw_baldauf_brdar",
    "test_lamb_wave": "pybella.tests.test_lamb_wave",
    "test_blending_warm_bubble": "pybella.tests.test_blending_warm_bubble",
    "test_blending_hydrostatic": "pybella.tests.test_blending_hydrostatic",
    "test_blending_swe": "pybella.tests.test_blending_swe",
    "test_unstable_lamb": "pybella.tests.test_unstable_lamb",
    "test_swe_vortex": "pybella.tests.test_swe_vortex",
    "test_sphere_swe_tc2": "pybella.tests.test_sphere_swe_tc2",
    "test_sphere_swe_tc2_global": "pybella.tests.test_sphere_swe_tc2_global",
    "test_sphere_swe_tc6": "pybella.tests.test_sphere_swe_tc6",
    "test_sphere_swe_tc6_global": "pybella.tests.test_sphere_swe_tc6_global",
    "test_sphere_gw": "pybella.tests.test_sphere_gw",
    "test_hj_baroclinic": "pybella.tests.test_hj_baroclinic",
    "test_hj_baroclinic_ridges": "pybella.tests.test_hj_baroclinic_ridges",
    "test_hj_baroclinic_global": "pybella.tests.test_hj_baroclinic_global",
    "test_hj_baroclinic_ridges_global": "pybella.tests.test_hj_baroclinic_ridges_global",
    "test_straka": "pybella.tests.test_straka",
    "test_agnesi_hydrostatic": "pybella.tests.test_agnesi_hydrostatic",
    "test_agnesi_3d": "pybella.tests.test_agnesi_3d",
    "test_schaer_ridge": "pybella.tests.test_schaer_ridge",
    "smoke_zvert": "pybella.tests.smoke_zvert",
    "smoke_agnesi": "pybella.tests.smoke_agnesi",
}
