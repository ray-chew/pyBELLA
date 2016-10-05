/*
 *  set_ghostnodes_p.h
 *  LowMach.¹
 *
 *  Created by WorkAccount on Sun Feb 22 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SET_GHOSTNODES_P_H
#define SET_GHOSTNODES_P_H

#include "kgrid.h"

void set_ghostnodes_p(
					  double* p, 
					  const NodeSpaceDiscr* node, 
					  const int ig);

#endif