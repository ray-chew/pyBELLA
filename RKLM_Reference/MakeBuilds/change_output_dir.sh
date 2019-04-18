#!/bin/bash

# Replacing output directory in target input file with current user's directory
#
# Usage (provided change_output_dir.sh is executable, otherwise chmod 755 change_output_dir.sh)
# source change_output_dir.sh test absolute_path
# where test is the (CamelCase) string contained in /Input/userdata_test.c and absolute_path is of the form /abs_path_to_output/
#
# Example:
# source change_output_dir.sh TravellingVortex_2D /home/benacchio/workspace/RKLM_Reference/
# 
# Note: to undo
# git checkout -- ../Input/userdata_test.c

DEF="OutputBaseFolder      = .*;"
DEF2="OutputBaseFolder      = \"$2\""

FILE_NAME="../Input/userdata_$1.c"

sed -i "s|$DEF|$DEF2;|g" $FILE_NAME
