/*
 *  set_ghostfaces_h_over_rho.h
 *  LowMach.¹
 *
 *  Created by WorkAccount on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef SET_GHOSTNODES_H_OVER_RHO
#define SET_GHOSTNODES_H_OVER_RHO

#include "kgrid.h"

void set_ghostfaces_h_over_rho(
							   double** hplus, 
							   double** hminus, 
							   double** lhplus, 
							   double** lhminus, 
							   const NodeSpaceDiscr* node, 
							   const int ig);

#endif