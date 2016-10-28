/*
 *  set_ghostcells_p2_S2.h
 *  LowMach.π
 *
 *  Created by WorkAccount on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef SET_GHOSTCELLS_P_H
#define SET_GHOSTCELLS_P_H

#include "Common.h"
#include "kgrid.h"
#include "variable.h"

void set_ghostcells_p2(
                       double* p,                        
                       const double* hplus[3],
					   const ElemSpaceDiscr* elem, 
					   const int ig);

#endif
