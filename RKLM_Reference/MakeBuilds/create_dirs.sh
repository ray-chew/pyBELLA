#!/bin/bash

# Creates directories for hdf output for test cases for paper. 
#
# Usage (provided it is executable, otherwise chmod 755 create_dirs.sh)
# source clean_dirs.sh

mkdir ../../hdf_output

for TEST_C in 'InternalWave_Baldauf' 'InternalWave_NH' 'InternalWave_H' 'InternalWave_P' 'Straka_3D' 'Straka_3D_400m' 'Straka_3D_200m' 'Straka_3D_100m' 'Straka_3D_50m' 'TravellingVortex_3D_48' 'TravellingVortex_3D_96' 'TravellingVortex_3D_192'
do
  mkdir ../../hdf_output/$TEST_C
  for VAR in 'advflux' 'dpdim' 'geopot' 'interface_data' 'p' 'qc' 'rho' 'rhoZp' 'rhs_nodes' 'second_projection_test' 'u' 'Z' 'buoy' 'drhoY' 'Graphics' 'lap_cells' 'p2_c' 'qr' 'rhoe' 'rhs' 'rhs_second' 'Srhs_cells' 'v' 'dp2_c' 'dY' 'hplus_nodes' 'lap_nodes' 'p2_nodes' 'qv' 'rhoY' 'rhs_cells' 'S' 'theta' 'w' 'dp2_nodes' 'fluxes' 'hydrostate' 'omega' 'p_prolongation' 'res' 'rhoZB' 'rhs_fg' 'S2_cells' 'time_series.txt' 'Y' 'psinc'
  do
    mkdir ../../hdf_output/$TEST_C/$VAR
  done
done