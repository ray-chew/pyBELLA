	/*
 *  variable_coefficient_poisson_nodes.c
 *  LowMach.π
 *
 *  Created by WorkAccount on Sat Feb 21 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include <float.h>

#include "error.h"
#include "SimpleUtilities.h"
#include "kgrid.h"
#include "mpv.h"
#include "thermodynamic.h"
#include "variable_coefficient_poisson_nodes.h"
#include "BiCGSTAB.h"
#include "laplacian_nodes.h"
#include "userdata.h"
#include "math_own.h"

#define DEBUG_OUTPUT 
#ifdef DEBUG_OUTPUT
#include "userdata.h"
#include "io.h"
extern User_Data ud;
#endif


/* ========================================================================== */

void set_periodic_data(double *p,
                       const NodeSpaceDiscr* node,
                       const int x_periodic,
                       const int y_periodic,
                       const int z_periodic
                       )
{
    const int igx = node->igx;
    const int igy = node->igy;
    const int igz = node->igz;
    const int icx = node->icx;
    const int icy = node->icy;
    const int icz = node->icz;
    
    int i, j, k, l, m, n_left, n_right;
    
    if (x_periodic) {
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                n_left  = m + igx;
                n_right = m + icx - igx - 1;
                p[n_right] = p[n_left];
            }
        }
    }
    if (node->ndim > 1 && y_periodic) {
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(i = igx; i < icx - igx; i++) {m = l + i;
                n_left  = m + igy * icx;
                n_right = m + (icy - igy - 1) * icx;
                p[n_right] = p[n_left];
            }
        }
    }
    if (node->ndim > 2 && z_periodic) {
        for(i = igx; i < icx - igx; i++) {l = i;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                n_left  = m + igz * icx * icy;
                n_right = m + (icz - igz - 1) * icx * icy;
                p[n_right] = p[n_left];
            }
        }
    }
}

/* ========================================================================== */

static double BiCGSTAB_MG_nodes(
								BiCGSTABData* data,
								const NodeSpaceDiscr* node,
								const ElemSpaceDiscr* elem,
								const double* hplus[3],
								const double* hcenter,
								const double* rhs,
								double* solution_io,
								const int x_periodic,
								const int y_periodic,
								const int z_periodic,
								const double dt) {
	
	/*
	 */
	
	double* r_0        = data->r_0;
	double* r_j        = data->r_j;
	double* p_j        = data->p_j;
	double* v_j        = data->v_j;
	double* s_j        = data->s_j;
	double* t_j        = data->t_j;
    double* rhs_prec   = data->help_vec;
	
	const int nc  = node->nc;
	const int igx = node->igx;
	const int igy = node->igy;
	const int igz = node->igz;
	const int icx = node->icx;
	const int icy = node->icy;
	const int icz = node->icz; 
	
	const double precision0  = data->precision;
	const int max_iterations = data->max_iterations;
	
    int cnt, i, j, k, l, m, n, cell_cnt;
	double rho1, alpha, beta, rho2, omega, sigma, tmp, tmp_local;
	double precision = precision0;
	
    const int imax = MAX_own(1,icx - igx - x_periodic);
    const int jmax = MAX_own(1,icy - igy - y_periodic);
    const int kmax = MAX_own(1,icz - igz - z_periodic);
    
	assert(data->size >= node->nc);                          /* Rupert: This could be dangerous. */
    
    
    double precon_inv_scale = precon_prepare(node, elem, hplus, hcenter, x_periodic, y_periodic, z_periodic);
    
    set_periodic_data(solution_io,	node, x_periodic, y_periodic, z_periodic);
    
    EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, solution_io, hplus, hcenter, x_periodic, y_periodic, z_periodic, v_j);

    precon_invert(rhs_prec, rhs, node, x_periodic, y_periodic, z_periodic);
    
    cell_cnt = 0;
	
	for(k = igz; k < kmax; k++) {l = k * icx * icy;
		for(j = igy; j < jmax; j++) {m = l + j * icx;
			for(i = igx; i < imax; i++) {n = m + i;
				r_j[n] = r_0[n] = rhs_prec[n] - v_j[n];
				p_j[n] = 0.0;
				v_j[n] = 0.0; 
                cell_cnt++;
			}
		}
	}
    

    tmp = 0.0;
    tmp_local = 0.0;
    for(k = igz; k < kmax; k++) {l = k * icx * icy;
        for(j = igy; j < jmax; j++) {m = l + j * icx;
            for(i = igx; i < imax; i++) {n = m + i;
                tmp += r_j[n] * r_j[n];
                tmp_local = MAX_own(tmp_local, fabs(r_j[n]));
            }
        }
    }

    alpha = omega = rho1 = 1.;
	tmp_local *= dt/(precon_inv_scale*precision);
    tmp = dt*sqrt(tmp/cell_cnt)/(precon_inv_scale*precision);
	
    printf(" iter = 0,  residual = %e,  local residual = %e,  gridsize = %d\n", tmp, tmp_local, nc);

	cnt = 0;
#ifdef DIV_CONTROL_LOCAL
    // while(tmp > 1.0 && tmp_local > 1.0 && cnt < 1 )
	while((tmp > 1.0 || tmp_local > 1.0) && cnt < max_iterations )
#else
    while(tmp > 1.0 && cnt < max_iterations)
#endif
    {
		rho2 = 0.0; 
		for(k = igz; k < kmax; k++) {l = k * icx * icy;
			for(j = igy; j < jmax; j++) {m = l + j * icx;
				for(i = igx; i < imax; i++) {n = m + i;
					rho2 += r_j[n] * r_0[n];
				}
			}
		}
		
		beta = (rho2 * alpha) / (rho1 * omega);
		
		for(k = igz; k < kmax; k++) {l = k * icx * icy;
			for(j = igy; j < jmax; j++) {m = l + j * icx;
				for(i = igx; i < imax; i++) {n = m + i;
					p_j[n] = r_j[n] + beta * (p_j[n] - omega * v_j[n]);
				}
			}
		}
		
        set_periodic_data(p_j,	node, x_periodic, y_periodic, z_periodic);
        
        EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, p_j, hplus, hcenter, x_periodic, y_periodic, z_periodic, v_j);

		sigma = 0.0; 
		for(k = igz; k < kmax; k++) {l = k * icx * icy;
			for(j = igy; j < jmax; j++) {m = l + j * icx;
				for(i = igx; i < imax; i++) {n = m + i;
					sigma += v_j[n] * r_0[n];
				}
			}
		}
		
		alpha = rho2 / sigma;
		
		for(k = igz; k < kmax; k++) {l = k * icx * icy;
			for(j = igy; j < jmax; j++) {m = l + j * icx;
				for(i = igx; i < imax; i++) {n = m + i;
					s_j[n] = r_j[n] - alpha * v_j[n];
				}
			}
		}
		
        set_periodic_data(s_j,	node, x_periodic, y_periodic, z_periodic);

        EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, s_j, hplus, hcenter, x_periodic, y_periodic, z_periodic, t_j);

		omega = 0.0; 
		tmp = 0.0; 
		for(k = igz; k < kmax; k++) {l = k * icx * icy;
			for(j = igy; j < jmax; j++) {m = l + j * icx;
				for(i = igx; i < imax; i++) {n = m + i;
					omega += s_j[n] * t_j[n];
					tmp += t_j[n] * t_j[n];
				}
			}
		}
		
		omega /= tmp;
		
		for(k = igz; k < kmax; k++) {l = k * icx * icy;
			for(j = igy; j < jmax; j++) {m = l + j * icx;
				for(i = igx; i < imax; i++) {n = m + i;
					solution_io[n] += alpha * p_j[n] + omega * s_j[n];
					r_j[n] = s_j[n] - omega * t_j[n];
				}
			}
		}
        
        tmp = 0.0;
        tmp_local = 0.0;
        for(k = igz; k < kmax; k++) {l = k * icx * icy;
            for(j = igy; j < jmax; j++) {m = l + j * icx;
                for(i = igx; i < imax; i++) {n = m + i;
                    tmp      += r_j[n] * r_j[n];
                    tmp_local = MAX_own(tmp_local, fabs(r_j[n]));
                }
            }
        }
        
		rho1 = rho2;
		tmp_local *= dt/(precon_inv_scale*precision);
        tmp = dt*sqrt(tmp/cell_cnt)/(precon_inv_scale*precision);
		cnt++;
		
		if(cnt % 100 == 0) printf(" iter = %d, residual = %e\n", cnt, tmp);  
        set_periodic_data(solution_io,	node, x_periodic, y_periodic, z_periodic);
	}
        
    // assert(cnt == 23);
	printf(" iter = %d,  residual = %e,  local residual = %e,  gridsize = %d\n", cnt, tmp, tmp_local, nc);  
	
	data->actual_iterations = cnt;
	
	return(tmp);
	
}

/* ========================================================================== */

void variable_coefficient_poisson_nodes(
                                        double *p2,
                                        const double *hplus[3],
                                        const double *hcenter,
                                        const double *rhs,
                                        const ElemSpaceDiscr* elem,
                                        const NodeSpaceDiscr* node,
                                        const int x_periodic,
                                        const int y_periodic,
                                        const int z_periodic,
                                        const double dt)
{
    extern User_Data ud;
    extern MPV* mpv;

    const int nc = node->nc;
    
    const double precision0 = ud.second_projection_precision;   /* CHECK: precision for residual in L2-Norm */
    const double local_precision0 = ud.second_projection_local_precision;   /* CHECK: precision for residual in L2-Norm */
    const int max_iter = ud.second_projection_max_iterations;  /* 6 for MG; max no. of BiCGSTAB-iterations in a smoothing step */
    const int outperiod = 50;                                  /* frequency of output from BiCGSTAB-iterator*/
    
    double tmp, precision = precision0, local_precision = local_precision0;
    
    BiCGSTABData* data;
    
    data = BiCGSTABData_new(nc, precision, local_precision, max_iter, outperiod);
    int maxit = data->max_iterations;
    data->max_iterations = ud.second_projection_max_iterations;
    tmp = BiCGSTAB_MG_nodes(data, node, elem, hplus, hcenter, rhs, p2, x_periodic, y_periodic, z_periodic, dt);
    printf("residual 2nd projection = %e * tol\n", tmp);
    
    data->max_iterations = maxit;
    
    BiCGSTABData_free(data); 
}
