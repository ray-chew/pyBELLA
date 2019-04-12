#!/bin/bash

# Replacing test case name in Makefiles with current test
#
# Usage (provided changetest.sh is executable, otherwise chmod 755 changetest.sh)
# source changetest.sh testname
# where the input file for testname is of the form Input/userdata_testname.c
#
# Example:
# source changetest.sh AcousticWave
# 
# Note: to undo
# git checkout -- Makefile ../Makefile Makefile.txt

for FILE in '../Makefile' 'Makefile' 'Makefile.txt'
do
sed -i "s/userdata_.*.c/userdata_$1.c/g" $FILE
sed -i "s/userdata_.*.o /userdata_$1.o /g" $FILE
done
