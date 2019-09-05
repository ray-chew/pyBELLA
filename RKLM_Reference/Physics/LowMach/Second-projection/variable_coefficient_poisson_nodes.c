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
        
    if (x_periodic) {
        for(int k = 0; k < icz; k++) {
            int l = k * icx * icy;
            for(int j = 0; j < icy; j++) {
                int m = l + j * icx;
                for (int ii=0; ii<=igx; ii++) {
                    int n_left  = m + igx + ii;
                    int n_right = m + icx - igx - 1 + ii;
                    p[n_right] = p[n_left];

                    n_left  = m + igx - ii;
                    n_right = m + icx - igx - 1 - ii;
                    p[n_left] = p[n_right];
                }
            }
        }
    }
    if (node->ndim > 1 && y_periodic) {
        for(int k = 0; k < icz; k++) {
            int l = k * icx * icy;
            for(int i = 0; i < icx; i++) {
                int m = l + i;
                for (int jj=0; jj<=igy; jj++) {
                    int n_left  = m + (igy + jj) * icx;
                    int n_right = m + (icy - igy - 1 + jj) * icx;
                    p[n_right] = p[n_left];

                    n_left  = m + (igy - jj) * icx;
                    n_right = m + (icy - igy - 1 - jj) * icx;
                    p[n_left] = p[n_right];
                }
            }
        }
    }
    if (node->ndim > 2 && z_periodic) {
        for(int i = 0; i < icx; i++) {
            int l = i;
            for(int j = 0; j < icy; j++) {
                int m = l + j * icx;
                for (int kk=0; kk<=igz; kk++) {
                    int n_left  = m + (igz + kk) * icx * icy;
                    int n_right = m + (icz - igz - 1 + kk) * icx * icy;
                    p[n_right] = p[n_left];

                    n_left  = m + (igz - kk) * icx * icy;
                    n_right = m + (icz - igz - 1 - kk) * icx * icy;
                    p[n_left] = p[n_right];
                }
            }
        }
    }
}

/* ========================================================================== */
BiCGSTABData* BiCGSTABData_new(
                               const int size,
                               const double precision,
                               const double local_precision,
                               const int max_iterations,
                               const int output_period) {
    
    BiCGSTABData* var = (BiCGSTABData*)malloc(sizeof(BiCGSTABData));
    
    var->size = size;
    
    var->r_0       = (double*)malloc(size * sizeof(double));
    var->r_j       = (double*)malloc(size * sizeof(double));
    var->p_j       = (double*)malloc(size * sizeof(double));
    var->v_j       = (double*)malloc(size * sizeof(double));
    var->s_j       = (double*)malloc(size * sizeof(double));
    var->t_j       = (double*)malloc(size * sizeof(double));
    var->help_vec  = (double*)malloc(size * sizeof(double));
    var->precision = precision;
    var->local_precision = local_precision;
    var->max_iterations = max_iterations;
    var->output_period = output_period;
    return var;
}

/* ========================================================================== */
void BiCGSTABData_free(BiCGSTABData* var) {
    free(var->r_0);
    free(var->r_j);
    free(var->p_j);
    free(var->v_j);
    free(var->s_j);
    free(var->t_j);
    free(var->help_vec);
    free(var);
}

/* ========================================================================== */

#if OUTPUT_LAP_NODES
static int lap_output_count = 0;
#define OUTPUT_RHS 1
#if OUTPUT_RHS
static int rhs_output_count = 0;
#endif
#endif

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
    // double precon_inv_scale = 1.0;
    
    set_periodic_data(solution_io, node, x_periodic, y_periodic, z_periodic);

    EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, solution_io, hplus, hcenter, x_periodic, y_periodic, z_periodic, v_j);

#if OUTPUT_LAP_NODES
    FILE *plapfile = NULL;
    char fn[120], fieldname[90];
    if (lap_output_count < 10) {
        sprintf(fn, "%s/lap_nodes/lap_nodes_00%d.hdf", ud.file_name, lap_output_count);
    } else if(lap_output_count < 100) {
        sprintf(fn, "%s/lap_nodes/lap_nodes_0%d.hdf", ud.file_name, lap_output_count);
    } else {
        sprintf(fn, "%s/lap_nodes/lap_nodes_%d.hdf", ud.file_name, lap_output_count);
    }
    sprintf(fieldname, "lap_nodes");    
    WriteHDF(plapfile, node->icx, node->icy, node->icz, node->ndim, v_j, fn, fieldname);
    lap_output_count++;
#endif

    precon_invert(rhs_prec, rhs, node, x_periodic, y_periodic, z_periodic);

#if OUTPUT_RHS
    FILE *prhsfile = NULL;
#ifndef OUTPUT_LAP_NODES
    char fn[120], fieldname[90];
#endif
    if (rhs_output_count < 10) {
        sprintf(fn, "%s/rhs_nodes/rhs_nodes_prec_00%d.hdf", ud.file_name, rhs_output_count);
    } else if(rhs_output_count < 100) {
        sprintf(fn, "%s/rhs_nodes/rhs_nodes_prec_0%d.hdf", ud.file_name, rhs_output_count);
    } else {
        sprintf(fn, "%s/rhs_nodes/rhs_nodes_prec_%d.hdf", ud.file_name, rhs_output_count);
    }
    sprintf(fieldname, "rhs_nodes_prec");    
    WriteHDF(prhsfile, node->icx, node->icy, node->icz, node->ndim, rhs_prec, fn, fieldname);
    rhs_output_count++;
#endif

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
	
    printf(" iter = 0, residual = %e, local residual = %e, gridsize = %d\n", tmp, tmp_local, nc);

	cnt = 0;
#if DIV_CONTROL_LOCAL
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
		
        set_periodic_data(p_j, node, x_periodic, y_periodic, z_periodic);

        EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, p_j, hplus, hcenter, x_periodic, y_periodic, z_periodic, v_j);

#if OUTPUT_LAP_NODES
        if (lap_output_count < 10) {
            sprintf(fn, "%s/lap_nodes/lap_nodes_00%d.hdf", ud.file_name, lap_output_count);
        } else if(lap_output_count < 100) {
            sprintf(fn, "%s/lap_nodes/lap_nodes_0%d.hdf", ud.file_name, lap_output_count);
        } else {
            sprintf(fn, "%s/lap_nodes/lap_nodes_%d.hdf", ud.file_name, lap_output_count);
        }
        sprintf(fieldname, "lap_nodes");    
        WriteHDF(plapfile, node->icx, node->icy, node->icz, node->ndim, v_j, fn, fieldname);
        lap_output_count++;
#endif

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
		
        set_periodic_data(s_j, node, x_periodic, y_periodic, z_periodic);
        EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, s_j, hplus, hcenter, x_periodic, y_periodic, z_periodic, t_j);

#if OUTPUT_LAP_NODES
        if (lap_output_count < 10) {
            sprintf(fn, "%s/lap_nodes/lap_nodes_00%d.hdf", ud.file_name, lap_output_count);
        } else if(lap_output_count < 100) {
            sprintf(fn, "%s/lap_nodes/lap_nodes_0%d.hdf", ud.file_name, lap_output_count);
        } else {
            sprintf(fn, "%s/lap_nodes/lap_nodes_%d.hdf", ud.file_name, lap_output_count);
        }
        sprintf(fieldname, "lap_nodes");    
        WriteHDF(plapfile, node->icx, node->icy, node->icz, node->ndim, t_j, fn, fieldname);
        lap_output_count++;
#endif

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
        set_periodic_data(solution_io, node, x_periodic, y_periodic, z_periodic);
	}
        
    // assert(cnt == 23);
	printf(" iter = %d, residual = %e, local residual = %e, gridsize = %d\n", cnt, tmp, tmp_local, nc);  
	
	data->actual_iterations = cnt;
	
	return(tmp);
	
}

/* ========================================================================== */
int shitty_count = 0;

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
    
    // double tmp_arr[53*53];
    // for (int ii = 0; ii < nc; ii ++) {
    //     tmp_arr[ii] = p2[ii];
    //     // printf("tmp_arr[ii] = %e", tmp_arr[ii]);
    // }

    // FILE *pnewfile = NULL;
    // char fn[120], fieldname[90];
    // sprintf(fn, "%s/pnew/pnew_00%d.hdf", ud.file_name, shitty_count);
    // sprintf(fieldname, "pnew");    
    // WriteHDF(pnewfile, node->icx, node->icy, node->icz, node->ndim, tmp_arr, fn, fieldname);
    // shitty_count += 1;

    BiCGSTABData_free(data); 
}
