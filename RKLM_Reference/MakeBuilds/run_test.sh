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

if [ -z "$1" ]
  then
    echo "Error: select test."
    exit 1;
fi

# Cleaning running output directories
source clean_hdf_output.sh
# Cleaning plotting output directories
source clean_hdf_output_test.sh $TESTNAME
# Changing test-case dependent instances in Makefiles
source changetest.sh $TESTNAME
# Changing output directory in userdata file to the one relative to current path
SRC_DIR=$(cd "$(dirname "$0")/../.."; pwd)
SRC_DIR+="/"
source change_output_dir.sh $TESTNAME $SRC_DIR


# Changing value of compressible flag according to modelstr, to be completed
#sed -i "s/ud->is_compressible   = 1/ud->is_compressible   = $MSV/g" "../Input/userdata_$TESTNAME.c"
#sed -i "s/ud->is_compressible   = 0/ud->is_compressible   = $MSV/g" "../Input/userdata_$TESTNAME.c"

# Cleaning the executables
make clean
# Compile the code
make

# Run the code 
./rklm 

# Copying results from running directory to directory used for plotting
cp -r ../../low_Mach_gravity_comp/* ../../hdf_output/$TESTNAME/
cp -r ../../low_Mach_gravity_psinc/* ../../hdf_output/$TESTNAME/

