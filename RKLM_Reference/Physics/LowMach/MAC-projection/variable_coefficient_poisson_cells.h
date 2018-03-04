/*
 *  variable_coefficient_poisson_cells.h
 *  LowMach.¹
 *
 *  Created by WorkAccount on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef VARIABLE_COEFFICIENT_POISSION_CELLS
#define VARIABLE_COEFFICIENT_POISSION_CELLS

#include "variable.h"

void variable_coefficient_poisson_cells(
										double *p2,
                                        double *rhs,
										const double *hplus[3],
										const double *hcenter,
										const ConsVars* Sol,
                                        const ElemSpaceDiscr* elem,
                                        const NodeSpaceDiscr* node);

#endif
