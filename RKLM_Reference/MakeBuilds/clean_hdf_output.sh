#!/bin/bash

# Cleans hdf output from comp and psinc folders 
#
# Usage (provided it is executable, otherwise chmod 755 clean_hdf_output.sh)
# source clean_hdf_output.sh

for FOLDER in 'comp' 'psinc'
do
  find ../../low_Mach_gravity_$FOLDER/*/ -name "*.hdf" -print0 | xargs -0 rm
done
