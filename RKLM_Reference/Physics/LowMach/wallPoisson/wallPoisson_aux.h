//
//  wallPoisson_aux.h
//  LowMach_MG_V
//
//  Created by Klein, Rupert on 29/11/14.
//  Copyright (c) 2014 FU Berlin. All rights reserved.
//

#ifndef __LowMach_MG_V__wallPoisson_aux__
#define __LowMach_MG_V__wallPoisson_aux__

#include <stdio.h>

double betaVal(double *xx);

void gradVal(double *xx, double* grad);

#endif /* defined(__LowMach_MG_V__wallPoisson_aux__) */
