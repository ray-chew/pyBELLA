/*
 *  set_ghostnodes_rhs.h
 *  LowMach.¹
 *
 *  Created by WorkAccount on Sun Feb 22 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SET_GHOSTNODES_RHS_H
#define SET_GHOSTNODES_RHS_H

#include "kgrid.h"

void set_ghostnodes_rhs(
						double* rhs, 
						const NodeSpaceDiscr* node, 
						const int ig);

#endif