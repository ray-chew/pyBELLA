#!/bin/bash

# Runs the suite of tests from RKLM low Mach fluid dynamics code
# for the paper "A semi-implicit numerical model for small-to-planetary scale atmospheric
# dynamics"

for TESTCASE in 'AcousticWave_low' 'AcousticWave_high' 'InternalWave_NH' 'InternalWave_H' 'InternalWave_P'
do
  ./run_test.sh $TESTCASE 
done