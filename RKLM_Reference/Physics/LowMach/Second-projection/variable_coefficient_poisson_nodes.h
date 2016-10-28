/*
 *  variable_coefficient_poisson.h
 *  LowMach.π
 *
 *  Created by WorkAccount on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#include "Common.h"

#ifndef VARIABLE_COEFFICIENT_POISSION_NODES
#define VARIABLE_COEFFICIENT_POISSION_NODES

void variable_coefficient_poisson_nodes(
										double *p2_argument,
										const double *hplus[3],
										const double *hcenter,
										const double *rhs,
										const int x_periodic,
										const int y_periodic, 
										const int z_periodic,
										const double dt);

#endif /* VARIABLE_COEFFICIENT_POISSION_NODES */ 
