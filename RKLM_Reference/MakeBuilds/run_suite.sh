#!/bin/bash

# Runs the suite of tests from RKLM low Mach fluid dynamics code
# for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
# dynamics"

for TESTCASE in 'TravellingVortex_3D' 'InternalWave_NH' 'InternalWave_H' 'InternalWave_P' 'Straka_3D_400m' 'Straka_3D_200m' 'Straka_3D_100m' 'Straka_3D_50m' 'InternalWave_Baldauf'
do
  ./run_test.sh $TESTCASE 
done