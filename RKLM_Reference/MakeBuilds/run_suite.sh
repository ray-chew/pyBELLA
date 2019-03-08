#!/bin/bash

# Runs the suite of tests from RKLM low Mach fluid dynamics code
# for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
# dynamics" and plots results
#

# Straka
#for TESTCASE in 'Straka_3D_400m' 'Straka_3D_200m' 'Straka_3D_100m' 'Straka_3D_50m'
# SK94
#for TESTCASE in 'InternalWave_P_psinc' 'InternalWave_P_hyd' 'InternalWave_P'
# Tests in paper
for TESTCASE in 'TravellingVortex_3D_768' 'TravellingVortex_3D_384' 'TravellingVortex_3D_192' 'TravellingVortex_3D_96' 'TravellingVortex_3D_48' #'InternalWave_NH' 'InternalWave_H' 'InternalWave_H_psinc' 'InternalWave_H_hyd' 'InternalWave_P' 'InternalWave_P_psinc' 'InternalWave_P_hyd' 'Straka_3D_400m' 'Straka_3D_200m' 'Straka_3D_100m' 'Straka_3D_50m' 'Straka_3D_25m'
#for TESTCASE in 'TravellingVortex_3D_1024' 'TravellingVortex_3D_512' 'TravellingVortex_3D_256' 'TravellingVortex_3D_128' 'TravellingVortex_3D_64'
do
  ./run_test.sh $TESTCASE 
done

# Plot results, make sure that tests in make plots are the ones executed above.
#cd ../../MatLab/ 
#matlab -nodesktop -r make_plots, quit
