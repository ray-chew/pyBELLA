//
//  wallPoisson_aux.c
//  LowMach_MG_V
//
//  Created by Klein, Rupert on 29/11/14.
//  Copyright (c) 2014 FU Berlin. All rights reserved.
//

#include "wallPoisson_aux.h"

/* ========================================================================== */

double betaVal(double *xx)
{
    /* coeff multiplying   grad p   in  div(coeff*grad p) = rhs  to be solved
     in the projection steps
     
     This is to be improved to allow a more flexible, user-defined, flow
     dependent assignment of coeff.
     */
    
    return 1;
}

/* ========================================================================== */

void gradVal(double *xx, double* grad)
{
    /* this is only relevant for cut cells */
    
     grad[0] = 1.0;
     grad[1] = 0.0;  
}
