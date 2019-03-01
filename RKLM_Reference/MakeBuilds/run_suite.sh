#!/bin/bash

# Runs the suite of tests from RKLM low Mach fluid dynamics code
# for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
# dynamics" and plots results
#

# Straka
#for TESTCASE in 'Straka_3D_400m' 'Straka_3D_200m' 'Straka_3D_100m' 'Straka_3D_50m'
# SK94
for TESTCASE in 'InternalWave_NH' 'InternalWave_H' 'InternalWave_P'
# Tests in paper
#for TESTCASE in 'TravellingVortex_3D_192' 'InternalWave_NH' 'InternalWave_H' 'InternalWave_P' 'Straka_3D_50m' 'InternalWave_Baldauf'
do
  ./run_test.sh $TESTCASE 
done

# Plot results, make sure that tests in make plots are the ones executed above.
cd ../../MatLab/ 
matlab -nodesktop -r make_plots, quit
