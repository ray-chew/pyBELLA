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

#include "ProjectionType.h"
#include "variable_coefficient_poisson_nodes.h"
#include "set_ghostnodes_p.h"
#include "set_ghostnodes_rhs.h"

#ifdef SOLVER_2_HYPRE_RUPE
#include "hypreInterface.h"
#endif

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
    if (y_periodic) {
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(i = igx; i < icx - igx; i++) {m = l + i;
                n_left  = m + igy * icx;
                n_right = m + (icx - igx - 1) * icx;
                p[n_right] = p[n_left];
            }
        }
    }
    if (z_periodic) {
        for(i = igx; i < icx - igx; i++) {l = i;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                n_left  = m + igz * icx * icy;
                n_right = m + (icz - igz - 1) * icx * icy;
                p[n_right] = p[n_left];
            }
        }
    }
}


#ifdef SOLVER_2_HYPRE
#ifdef SOLVER_2_HYPRE_RUPE

#include "HYPRE_sstruct_ls.h"
#include "hypreInterface.h"

/* ========================================================================== */

void matrixInnerElements(const double*          hplus[3],
                         const double*          hcenter,
                         const double*          hgrav,
                         const int              x_periodic,
                         const int              y_periodic,
                         const int              z_periodic,
                         HYPRE_SStructMatrix    A,
                         HYPRE_SStructVector    b)
{    

    extern User_Data ud;
    extern MPV* mpv;
    const ElemSpaceDiscr* elem = mpv->Level[0]->elem;
    const NodeSpaceDiscr* node = mpv->Level[0]->node;

    const int icxe = elem->icx;
    const int icye = elem->icy;
    const int igxe = elem->igx;
    const int igye = elem->igy;
        
    const double dV   = 1.0;                         /* node->dx*node->dy; */
    const double dxdy = node->dx / node->dy / dV;
    const double dydx = node->dy / node->dx / dV;
    
    const double nine_pt  = (0.25 * (1.0 + P2_DIAGONAL_FIVE_POINT)) * P2_FULL_STENCIL;
    const double wcloseby = 0.5 * (1.0-nine_pt);
    const double wfaraway = 0.5 * (0.0+nine_pt);
    
    const int stridexn = (ud.bdrytype_min[0] == PERIODIC ? elem->icx : node->icx);
    
    for (int j = igye; j < icye-igye; j++) {
        
        int me = j*icxe;
        int jn = j-igye; 
        
        for (int i = igxe; i < icxe-igxe; i++) {
            
            int ne   = me+i;
            int in   = i-igxe;
            int inp1 = (in+1)%stridexn; 
                        
            double IntDpDn[4][4];
            
            double M[4][4] = {  {-wcloseby, wcloseby, wfaraway,-wfaraway},
                                {-wfaraway, wfaraway, wcloseby,-wcloseby},
                                {-wcloseby,-wfaraway, wfaraway, wcloseby},
                                {-wfaraway,-wcloseby, wcloseby, wfaraway}};
            
            /* for a single primary cell obtain weights for its nodal
               pressures for determining the following edge integrals
               on boundaries of dual cells:
             
             
                (i,j+1)---------(i+1,j+1)
                      |         |
                      |         |
                      |         |
                      |         |
                  (i,j)---------(i+1,j)

             
             
                -----   -----   -----   -----
               |     | |     | |  |  | |  |  |
               |--   | |   --| |   --| |--   |
               |  |  | |  |  | |     | |     |
                -----   -----   -----   -----
                IDP0    IDP1    IDP2    IDP3
               
               where IDPi translates to IntDpDn[i][]
             
             */
            for (int k = 0; k < 4; k++)
            {
                M[0][k] *= dydx * hplus[0][ne];
                M[1][k] *= dydx * hplus[0][ne];
                M[2][k] *= dxdy * hplus[1][ne];
                M[3][k] *= dxdy * hplus[1][ne];
                                                      // target nodes
                IntDpDn[0][k] = + M[0][k] + M[2][k];  // (i,j)
                IntDpDn[1][k] = - M[0][k] + M[3][k];  // (i+1,j)
                IntDpDn[2][k] = - M[1][k] - M[3][k];  // (i+1,j+1)
                IntDpDn[3][k] = + M[1][k] - M[2][k];  // (i,j+1)
            }
                        
            int entries[4][4] = {{0, 2, 8, 4},
                                 {1, 0, 4, 7},
                                 {5, 3, 0, 1},
                                 {3, 6, 2, 0}};
            
            int ij[4][2] = {{in,jn}, {inp1,jn}, {inp1,jn+1}, {in,jn+1}};
            
            for (int k = 0; k < 4; k++)
            {
                HYPRE_SStructMatrixAddToValues(A, 0, ij[k], 0, 4, entries[k], IntDpDn[k]);
            }
        }
    }
}

/* ========================================================================== */

void map_rhs_to_HYPRE_format(HYPRE_SStructVector b, 
                             const double *rhs, 
                             const NodeSpaceDiscr *node){

    const int icxn = node->icx;
    const int icyn = node->icy;
    const int igxn = node->igx;
    const int igyn = node->igy;
        
    for (int j=igyn; j<icyn-igyn; j++) {
        int mn       = j*icxn;
        for (int i=igxn; i<icxn-igxn; i++) {
            int nn       = mn+i;
            int ij[] = {i-igxn , j-igyn};
            HYPRE_SStructVectorSetValues(b, 0, ij, 0, &rhs[nn]);
        }
    }
}

/* ========================================================================== */

void variable_coefficient_poisson_nodes(double *p2,
                                        const double *hplus[3],
                                        const double *hcenter,
                                        const double *hgrav,
                                        const double *rhs,
                                        const int x_periodic,
                                        const int y_periodic,
                                        const int z_periodic,
                                        const double dt)
{
    extern User_Data ud;
    extern MPV* mpv;
    const NodeSpaceDiscr* node = mpv->Level[0]->node;

    const int icx = node->icx;
    const int icy = node->icy;
    const int igx = node->igx;
    const int igy = node->igy;
    
    const int nx  = icx-2*igx;
    const int ny  = icy-2*igy;
    
    HYPRE_SStructMatrix A;
    HYPRE_SStructVector b;
    HYPRE_SStructVector x;
    
    initialize_Hypre_Axb(nx, ny, WRONG, &A, &x, &b);
    
    /* build matrix */
    matrixInnerElements(hplus, hcenter, hgrav, x_periodic, y_periodic, z_periodic, A, b);

    /* map rhs to HYPRE-format */
    map_rhs_to_HYPRE_format(b, rhs, node);
    
    /* solve */
    solve_with_Hypre(A, x, b);
    
    HYPRE_SStructVectorGather(x);
    
#ifdef HYPRE_PRINT
    HYPRE_SStructVectorPrint("x", x, 0);
#endif

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            int ij[] = {i , j};
            int nn   = (j+igy)*icx + (i+igx);
            HYPRE_SStructVectorGetValues(x, 0, ij, 0, &p2[nn]);
        }
    }
    
    HYPRE_SStructVectorDestroy(x);
    HYPRE_SStructVectorDestroy(b);
    HYPRE_SStructMatrixDestroy(A);
}

#endif /* SOLVER_2_HYPRE_RUPE */
#else  /* SOLVER_2_HYPRE */

/* ========================================================================== */

static double BiCGSTAB_MG_nodes(
								BiCGSTABData* data,
								const NodeSpaceDiscr* node,
								const ElemSpaceDiscr* elem,
								const double* hplus[3],
								const double* hcenter,
								const double* hgrav,
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
    double* r_j_unprec = data->help_vec2;
	
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
	
	assert(data->size >= node->nc);                          /* Rupert: This could be dangerous. */
    
    
    precon_prepare(node, elem, hplus, hcenter, hgrav, x_periodic, y_periodic, z_periodic);
    
    set_periodic_data(solution_io,	node, x_periodic, y_periodic, z_periodic);
    
    EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, solution_io, hplus, hcenter, hgrav, x_periodic, y_periodic, z_periodic, v_j);

    precon_invert(rhs_prec, rhs, node);
    
    cell_cnt = 0;
	
	for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icx * icy;
		for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx;
			for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
				r_j[n] = r_0[n] = rhs_prec[n] - v_j[n];
				p_j[n] = 0.0;
				v_j[n] = 0.0; 
                cell_cnt++;
			}
		}
	}
    
    precon_apply(r_j_unprec, r_j, node);

#if 0 /* TEST */
    {
        double err      = 0.0;
        double r_j_norm = 0.0;
        double *r_j_prec = (double*)malloc(node->nc*sizeof(double));
        memset(r_j_prec, 0.0, sizeof(double));
        precon_invert(r_j_prec, r_j_unprec, node);
        for (int nn=0; nn<node->nc; nn++) {
            r_j_norm += r_j[nn]*r_j[nn];
            err += (r_j_prec[nn]-r_j[nn])*(r_j_prec[nn]-r_j[nn]);
        }
        err      = sqrt(err/node->nc);
        r_j_norm = sqrt(r_j_norm/node->nc);
        printf("err, r_j_norm = %e, %e\n", err, r_j_norm);
    }
#endif
    
    tmp = 0.0;
    tmp_local = 0.0;
    for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icx * icy;
        for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
                tmp += r_j_unprec[n] * r_j_unprec[n];
                tmp_local = MAX_own(tmp_local, fabs(r_j_unprec[n]));
            }
        }
    }

    alpha = omega = rho1 = 1.;
	tmp_local *= 0.5*dt/precision;
    tmp = 0.5*dt*sqrt(tmp/cell_cnt)/precision;
	
    printf(" iter = 0,  residual = %e,  local residual = %e,  gridsize = %d\n", tmp, tmp_local, nc);

	cnt = 0;
#ifdef DIV_CONTROL_LOCAL
	while(tmp > 1.0 && tmp_local > 1.0 && cnt < max_iterations)
    {
#else
    while(tmp > 1.0 && cnt < max_iterations)
	{
#endif
		rho2 = 0.0; 
		for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icx * icy;
			for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx;
				for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
					rho2 += r_j[n] * r_0[n];
				}
			}
		}
		
		beta = (rho2 * alpha) / (rho1 * omega);
		
		for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icx * icy;
			for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx;
				for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
					p_j[n] = r_j[n] + beta * (p_j[n] - omega * v_j[n]);
				}
			}
		}
		
        set_periodic_data(p_j,	node, x_periodic, y_periodic, z_periodic);
        
        EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, p_j, hplus, hcenter, hgrav, x_periodic, y_periodic, z_periodic, v_j);

		sigma = 0.0; 
		for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icx * icy;
			for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx;
				for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
					sigma += v_j[n] * r_0[n];
				}
			}
		}
		
		alpha = rho2 / sigma;
		
		for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icx * icy;
			for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx;
				for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
					s_j[n] = r_j[n] - alpha * v_j[n];
				}
			}
		}
		
        set_periodic_data(s_j,	node, x_periodic, y_periodic, z_periodic);

        EnthalpyWeightedLap_Node_bilinear_p_scatter(node, elem, s_j, hplus, hcenter, hgrav, x_periodic, y_periodic, z_periodic, t_j);

		omega = 0.0; 
		tmp = 0.0; 
		for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icx * icy;
			for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx;
				for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
					omega += s_j[n] * t_j[n];
					tmp += t_j[n] * t_j[n];
				}
			}
		}
		
		omega /= tmp;
		
		for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icx * icy;
			for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx;
				for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
					solution_io[n] += alpha * p_j[n] + omega * s_j[n];
					r_j[n] = s_j[n] - omega * t_j[n];
				}
			}
		}

        precon_apply(r_j_unprec, r_j, node);
        
        tmp = 0.0;
        tmp_local = 0.0;
        for(k = igz; k < icz - igz - z_periodic; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy - y_periodic; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx - x_periodic; i++) {n = m + i;
                    tmp      += r_j_unprec[n] * r_j_unprec[n];
                    tmp_local = MAX_own(tmp_local, fabs(r_j_unprec[n]));
                }
            }
        }

		rho1 = rho2;
		tmp_local *= 0.5*dt/precision;
        tmp = 0.5*dt*sqrt(tmp/cell_cnt)/precision;
		cnt++;
		
		if(cnt % 100 == 0) printf(" iter = %d, residual = %e\n", cnt, tmp);  
        set_periodic_data(solution_io,	node, x_periodic, y_periodic, z_periodic);
	}
	printf(" iter = %d,  residual = %e,  local residual = %e,  gridsize = %d\n", cnt, tmp, tmp_local, nc);  
	
	data->actual_iterations = cnt;
	
	return(tmp);
	
}

/* ========================================================================== */

void variable_coefficient_poisson_nodes(
                                        double *p2,
                                        const double *hplus[3],
                                        const double *hcenter,
                                        const double *hgrav,
                                        const double *rhs,
                                        const int x_periodic,
                                        const int y_periodic,
                                        const int z_periodic,
                                        const double dt)
{
    extern User_Data ud;
    extern MPV* mpv;
    const NodeSpaceDiscr* node = mpv->Level[0]->node;
    const ElemSpaceDiscr* elem = mpv->Level[0]->elem;
    
    const int nc = mpv->Level[0]->node->nc;
    
    const double precision0 = ud.second_projection_precision;   /* CHECK: precision for residual in L2-Norm */
    const double local_precision0 = ud.second_projection_local_precision;   /* CHECK: precision for residual in L2-Norm */
    const int max_iter = ud.second_projection_max_iterations;  /* 6 for MG; max no. of BiCGSTAB-iterations in a smoothing step */
    const int outperiod = 50;                                  /* frequency of output from BiCGSTAB-iterator*/
    
    double tmp, precision = precision0, local_precision = local_precision0;
    
    BiCGSTABData* data;
    
    data = BiCGSTABData_new(nc, precision, local_precision, max_iter, outperiod);
    int maxit = data->max_iterations;
    data->max_iterations = ud.second_projection_max_iterations;
    tmp = BiCGSTAB_MG_nodes(data, node, elem, hplus, hcenter, hgrav, rhs, p2, x_periodic, y_periodic, z_periodic, dt);
    printf("residual 2nd projection = %e * tol\n", tmp);
    
    data->max_iterations = maxit;
    
    BiCGSTABData_free(data); 
}

#endif /* SOLVER_2_HYPRE */
