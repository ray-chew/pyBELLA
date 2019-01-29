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


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
typedef struct {
    
    int size;
    
    double* r_0;
    double* r_j;
    double* p_j;
    double* v_j;
    double* s_j;
    double* t_j;
    double* help_vec;
    
#ifdef SOLVER_CR2
    double* Lr_0;
#endif /* SOLVER_CR2 */
    
    double precision;
    double local_precision;
    int max_iterations;
    int actual_iterations;
    int output_period;
    
} BiCGSTABData;


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
BiCGSTABData* BiCGSTABData_new(
                               const int size,
                               const double precision,
                               const double local_precision,
                               const int max_iterations,
                               const int output_period);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void BiCGSTABData_free(BiCGSTABData* var);

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void variable_coefficient_poisson_nodes(
										double *p2_argument,
										const double *hplus[3],
										const double *hcenter,
										const double *rhs,
                                        const ElemSpaceDiscr* elem,
                                        const NodeSpaceDiscr* node,
										const int x_periodic,
										const int y_periodic, 
										const int z_periodic,
										const double dt);

#endif /* VARIABLE_COEFFICIENT_POISSION_NODES */ 
