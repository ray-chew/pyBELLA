# Initial-condition registry: -ic keys -> module paths providing UserData + sol_init.
#
# Legacy pre-restructure ICs (acoustic waves, SWE variants, baroclinic instability, ...)
# were removed when the old RKLM_Python/inputs tree was dropped; they remain recoverable
# from the git tag `archive/full_coriolis` and are mostly small tweaks of the cases below.
IC_MODULES = {
    "rb": "pybella.inputs.rising_bubble",
    "test_travelling_vortex": "pybella.tests.test_travelling_vortex",
    "test_travelling_vortex_3d_coriolis": "pybella.tests.test_travelling_vortex_3d_coriolis",
    "test_internal_long_wave": "pybella.tests.test_internal_long_wave",
    "test_igw_baldauf_brdar": "pybella.tests.test_igw_baldauf_brdar",
    "test_lamb_wave": "pybella.tests.test_lamb_wave",
    "test_blending_warm_bubble": "pybella.tests.test_blending_warm_bubble",
    "test_unstable_lamb": "pybella.tests.test_unstable_lamb",
    "test_swe_vortex": "pybella.tests.test_swe_vortex",
    "test_straka": "pybella.tests.test_straka",
    "smoke_zvert": "pybella.tests.smoke_zvert",
    "smoke_agnesi": "pybella.tests.smoke_agnesi",
}
