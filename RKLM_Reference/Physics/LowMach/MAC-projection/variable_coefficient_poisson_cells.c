/*
 *  variable_coefficient_poisson_cells.c
 *  LowMach.π
 *
 *  Created by WorkAccount on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <float.h>
#include <string.h>

#include "Common.h"

#include "error.h"


#include "variable_coefficient_poisson_cells.h"
#include "mpv.h"
#include "BiCGSTAB.h"
#include "laplacian_cells.h"
#include "userdata.h"
#include "thermodynamic.h"
#include "io.h"
#include "math_own.h"
#include "variable.h"
#include "kgrid.h"

/* ========================================================================== */

void variable_coefficient_poisson_cells(
										double *p2,
                                        double *rhs,
										const double *hplus[3],
                                        const double *hcenter,
										const ConsVars* Sol,
                                        const ElemSpaceDiscr* elem,
                                        const NodeSpaceDiscr* node)
{
	extern User_Data ud;
	extern MPV* mpv;
			
	const int nc = elem->nc;
	
	const double precision = ud.flux_correction_precision;  
	const double local_precision = ud.flux_correction_local_precision; 
	const int max_iter = ud.flux_correction_max_iterations; 
	const int outperiod = 50;                               
	
	double tmp;
	
	BiCGSTABData* data;

    memset(p2, 0.0, elem->nc*sizeof(double));

	data = BiCGSTABData_new(nc, precision, local_precision, max_iter, outperiod);
		
	tmp = SOLVER(data, elem, node, hplus, hcenter, Sol, mpv, mpv->dt, rhs, p2);
    
    printf("residual 1st projection = %e * tol\n", tmp);
	    
	BiCGSTABData_free(data); 
}

