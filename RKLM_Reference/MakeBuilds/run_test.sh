#!/bin/bash

# Run test
#
# Cleans up directory, changes files to work with test, builds and runs the testname test.  
# Usage (provided changetest.sh is executable, otherwise chmod 755 changetest.sh)
#
# source run_test.sh testname
#
# testname, select one of the following:
#           'Travelling-Vortex'   : see Kadioglu et al. 2008
#           'AcousticWave_low'    : see Vater 2013, Benacchio 2014
#           'AcousticWave_high'   : see Vater 2013, Benacchio 2014
#           'Straka_3D'           : density current, see Straka et al. 1993
#           'InternalWave_NH'     : nonhydrostatic inertia-gravity wave, see
#                                   Skamarock-Klemp 1994
#           'InternalWave_H'      : hydrostatic inertia-gravity wave, see
#                                   Skamarock-Klemp 1994
#           'InternalWave_P'      : planetary inertia-gravity wave, new test
#           'InternalWave_Baldauf', inertia-gravity wave test with rotation, see Baldauf-Brdar 2013
#
#
#

TESTNAME=$1

if [-z "$1"]
then
  echo "Error: select test."
  exit 1;
fi

source clean_hdf_output.sh
source clean_hdf_output_test.sh $TESTNAME
source changetest.sh $TESTNAME

# Changing value of compressible flag according to modelstr, to be completed
#sed -i "s/ud->is_compressible   = 1/ud->is_compressible   = $MSV/g" "../Input/userdata_$TESTNAME.c"
#sed -i "s/ud->is_compressible   = 0/ud->is_compressible   = $MSV/g" "../Input/userdata_$TESTNAME.c"

make clean
make
./rklm 

cp -r ../../low_Mach_gravity_comp/* ../../hdf_output/$TESTNAME/
cp -r ../../low_Mach_gravity_psinc/* ../../hdf_output/$TESTNAME/