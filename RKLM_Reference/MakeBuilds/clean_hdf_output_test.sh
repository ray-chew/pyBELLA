#!/bin/bash

# Cleans hdf output from output folders for specific test case
#
# Usage (provided it is executable, otherwise chmod 755 clean_hdf_output_test.sh)
# source clean_hdf_output_test.sh testname
#

TESTNAME=$1

find ../../hdf_output/$TESTNAME/*/ -name "*.hdf" -print0 | xargs -0 rm

