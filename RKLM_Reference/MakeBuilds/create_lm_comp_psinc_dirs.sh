#!/bin/bash

# Creates directories for hdf output for RKLM code. 
#
# Usage (provided it is executable, otherwise chmod 755 create_lm_comp_psinc_dirs.sh)
# source clean_lm_comp_psinc_dirs.sh

mkdir ../../low_Mach_gravity_comp
mkdir ../../low_Mach_gravity_psinc

for VAR in 'advflux' 'dpdim' 'geopot' 'interface_data' 'p' 'qc' 'rho' 'rhoZp' 'rhs_nodes' 'second_projection_test' 'u' 'Z' 'buoy' 'drhoY' 'Graphics' 'lap_cells' 'p2_c' 'qr' 'rhoe' 'rhs' 'rhs_second' 'Srhs_cells' 'v' 'dp2_c' 'dY' 'hplus_nodes' 'lap_nodes' 'p2_nodes' 'qv' 'rhoY' 'rhs_cells' 'S' 'theta' 'w' 'dp2_nodes' 'fluxes' 'hydrostate' 'omega' 'p_prolongation' 'res' 'rhoZB' 'rhs_fg' 'S2_cells' 'time_series.txt' 'Y' 'psinc'
do
  mkdir ../../low_Mach_gravity_comp/$VAR
  mkdir ../../low_Mach_gravity_psinc/$VAR
done
