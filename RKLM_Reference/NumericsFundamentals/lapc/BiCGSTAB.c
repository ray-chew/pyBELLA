/*******************************************************************************
 File:   BiCGSTAB.c
 Author: Nicola
 Date:   Sat Mar  7 19:26:55 CET 1998
 
 Notice: C version of Andreas Meister's C++ BiCGSTAB class
 *******************************************************************************/
#include "BiCGSTAB.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "Common.h"
#include "BiCGSTAB.h"
#include "kgrid.h"
#include "variable.h"
#include "mpv.h"
#include "math_own.h"
#include "laplacian_cells.h"
#include "boundary.h"


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
    
#ifdef SOLVER_CR2
    var->Lr_0      = (double*)malloc(size * sizeof(double));
#endif /* SOLVER_CR2 */
    
    var->precision = precision;
    var->local_precision = local_precision;
    var->max_iterations = max_iterations;
    var->output_period = output_period;
    
    return var;
}

void BiCGSTABData_free(BiCGSTABData* var) {
    free(var->r_0);
    free(var->r_j);
    free(var->p_j);
    free(var->v_j);
    free(var->s_j);
    free(var->t_j);
    free(var->help_vec);
    
#ifdef SOLVER_CR2
    free(var->Lr_0);
#endif /* SOLVER_CR2 */
    
    free(var);
}

#ifdef SOLVER_1_CR2
/* ========================================================================== */
/* Piotr's Conjugate Residual scheme */

double SOLVER(
              BiCGSTABData* data,
              const ElemSpaceDiscr* elem,
              const NodeSpaceDiscr* node,
              const double* hplus[3],
              const double* hcenter,
              const ConsVars* Sol,
              const MPV* mpv,
              const double dt,
              double* rhs,
              double* p2) {
    
    /*
     */
    
    double* r_m = data->r_0;
    double* r_0 = data->r_j;
    double* r_p = data->p_j;
    double* phi_m = data->v_j;
    double* phi_0 = data->s_j;
    double* phi_p = data->t_j;
    double* Lr_0  = data->help_vec;
    double* pshuffle;
    
    const int nc  = elem->nc;
    const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
    
    const double precision = data->precision;
    const double local_precision = data->local_precision;
    const int max_iterations = data->max_iterations;
    
    int cnt, cell_cnt, i, j, k, l, m, n;
    double AA, AB, AC, BB, BC, tmp, tmp_local, alpha, beta;
    
    assert(data->size >= elem->nc);                          /* Rupert: This could be dangerous. */
    
    double precon_inv_scale = precon_c_prepare(node, elem, hplus, hcenter);
    
    /* initialization */
    /* initial residual; intermediate abuse of the field */
    set_ghostcells_p2(p2, elem, 1);
    EnthalpyWeightedLap_bilinear_p(elem, node, p2, hplus, hcenter, Sol, mpv, dt, r_0);
    
    precon_c_invert(rhs, rhs, elem);
    
    for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                r_0[n]  -= rhs[n];
                phi_0[n] = p2[n];
                r_m[n]   = phi_m[n] = 0.0;
            }
        }
    }
    
    /* initialization of iteration coefficients; residual norm */
    set_ghostcells_p2(r_0, elem, 1);
    EnthalpyWeightedLap_bilinear_p(elem, node, r_0, hplus, hcenter, Sol, mpv, dt, Lr_0);
    
    /* norm of residuum */
    tmp = 0.0;
    tmp_local = 0.0;
    cell_cnt = 0;
    AA = 0.0;
    BB = 0.0;
    for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                tmp_local = MAX_own(tmp_local, fabs(r_0[n]));
                tmp += r_0[n] * r_0[n];
                AA  += r_0[n]  * Lr_0[n];
                BB  += Lr_0[n] * Lr_0[n];
                cell_cnt++;
            }
        }
    }
    tmp_local *= 0.5*dt/(precon_inv_scale*local_precision);
    tmp = 0.5*dt*sqrt(tmp/cell_cnt)/(precon_inv_scale*precision);
    
    alpha = 1.0;
    beta  = - AA / BB;
    
    /* iteration loop */
    cnt = 0;
    printf(" iter = 0,  residual = %e,  local residual = %e,  gridsize = %d\n", tmp, tmp_local, nc);
    
    /* while(tmp > 1.0 && tmp_local > 1.0 && cnt < max_iterations) */
    while(tmp > 1.0 && cnt < max_iterations)
    {
        
        /* iteration step */
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    phi_p[n] = alpha * (phi_0[n] - phi_m[n]) + phi_m[n] + beta * r_0[n];
                    r_p[n]   = alpha * (r_0[n]   - r_m[n])   + r_m[n]   + beta * Lr_0[n];
                }
            }
        }
        
        /* norm of residuum */
        tmp_local = 0.0;
        tmp       = 0.0;
        cell_cnt  = 0;
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    tmp_local = MAX_own(tmp_local, fabs(r_p[n]));
                    tmp      += (r_p[n]) * (r_p[n]);
                    cell_cnt++;
                }
            }
        }
        tmp_local *= 0.5*dt/(precon_inv_scale*local_precision);
        tmp = 0.5*dt*sqrt(tmp/cell_cnt)/(precon_inv_scale*precision);
        
        /* finish loop and reshuffle pointers */
        pshuffle = r_m;
        r_m = r_0;
        r_0 = r_p;
        r_p = pshuffle;
        
        pshuffle = phi_m;
        phi_m = phi_0;
        phi_0 = phi_p;
        phi_p = pshuffle;
        
        /* AA */
        AA = 0.0;
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    AA += (r_0[n] - r_m[n]) * (r_0[n] - r_m[n]);
                }
            }
        }
        
        /* AB */
        set_ghostcells_p2(r_0, elem, 1);
        EnthalpyWeightedLap_bilinear_p(elem, node, r_0, hplus, hcenter, Sol, mpv, dt, Lr_0);
        
        AB = 0.0;
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    AB += (r_0[n] - r_m[n]) * Lr_0[n];
                }
            }
        }
        
        /* AC */
        AC = 0.0;
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    AC += (r_0[n] - r_m[n]) * r_m[n];
                }
            }
        }
        
        /* BB  */
        BB = 0.0;
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    BB  += Lr_0[n] * Lr_0[n];
                }
            }
        }
        
        /* BC */
        BC = 0.0;
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    BC += Lr_0[n] * r_m[n];
                }
            }
        }
        
        /* recursion coefficients (alpha == Piotr's gamma) */
        alpha = (BC*AB - AC*BB) / (AA*BB - AB*AB);
        beta  = (AB*AC - AA*BC) / (AA*BB - AB*AB);
        
        cnt++;
        if(cnt % 100 == 0) {
            printf(" iter = %d,  residual = %e,  local residual = %e\n", cnt, tmp, tmp_local);
        }
    }
    
    /* store result into i/o-field */
    for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                p2[n] = phi_0[n];
            }
        }
    }
    
    set_ghostcells_p2(p2, elem, 1);
    
    printf(" iter = %d,  residual = %e,  local residual = %e,  gridsize = %d\n", cnt, tmp, tmp_local, nc);
    
    data->actual_iterations = cnt;
    
    return(MIN_own(tmp,tmp_local));
    
}

#else /* SOLVER_1_CR2 */

/* ========================================================================== */
/* BiCGSTAB */

double SOLVER(
              BiCGSTABData* data,
              const ElemSpaceDiscr* elem,
              const NodeSpaceDiscr* node,
              const double* hplus[3],
              const double* hcenter,
              const ConsVars* Sol,
              const MPV* mpv,
              const double dt,
              double* rhs,
              double* p2) {
    
    double* r_0 = data->r_0;
    double* r_j = data->r_j;
    double* p_j = data->p_j;
    double* v_j = data->v_j;
    double* s_j = data->s_j;
    double* t_j = data->t_j;
    double* rhs_prec   = data->help_vec;
    
    const int nc  = elem->nc;
    const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
    
    const double precision = data->precision;
    const double local_precision = data->local_precision;
    const int max_iterations = data->max_iterations;
    
    int cnt, i, j, k, l, m, n, cell_cnt;
    double rho1, alpha, beta, rho2, omega, sigma, tmp, tmp_local;
    
    assert(data->size >= elem->nc);                          /* Rupert: This could be dangerous. */
    
    double precon_inv_scale = precon_c_prepare(node, elem, hplus, hcenter);
    
    set_ghostcells_p2(p2, elem, 1);
    EnthalpyWeightedLap_bilinear_p(elem, node, p2, hplus, hcenter, Sol, mpv, dt, v_j);
    
    precon_c_invert(rhs_prec, rhs, elem);
    
    tmp_local = 0.0;
    cell_cnt = 0;
    
    for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                r_j[n]  = r_0[n]  = rhs_prec[n]   - v_j[n];
                p_j[n]  = 0.0;
                v_j[n]  = 0.0;
                cell_cnt++;
            }
        }
    }
    
    tmp = 0.0;
    tmp_local = 0.0;
    for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                tmp += r_j[n] * r_j[n];
                tmp_local = MAX_own(tmp_local, fabs(r_j[n]));
            }
        }
    }
    
    alpha = omega = rho1 = 1.;
    tmp_local *= 0.5*dt/(precon_inv_scale*local_precision);
    tmp = 0.5*dt*sqrt(tmp/cell_cnt)/(precon_inv_scale*precision);
    
    cnt = 0;
    printf(" iter = 0,  residual = %e,  local residual = %e,  gridsize = %d\n", tmp, tmp_local, nc);
    
    /* DIV_CONTROL_LOCAL  */
    while(tmp > 1.0 && tmp_local > 1.0 && cnt < max_iterations)
    {
        /*
         while(tmp > 1.0 && cnt < max_iterations)
         {
         */
        rho2 = 0.0;
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    rho2 += r_j[n] * r_0[n];
                }
            }
        }
        
        beta = (rho2 * alpha) / (rho1 * omega);
        
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    p_j[n]  = r_j[n]  + beta * (p_j[n]  - omega * v_j[n]);
                }
            }
        }
        
        set_ghostcells_p2(p_j, elem, 1);
        EnthalpyWeightedLap_bilinear_p(elem, node, p_j, hplus, hcenter, Sol, mpv, dt, v_j);
        
        sigma = 0.0; 
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    sigma += v_j[n] * r_0[n];
                }
            }
        }
        
        alpha = rho2 / sigma;
        
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    s_j[n]  = r_j[n]  - alpha * v_j[n];
                }
            }
        }
        
        set_ghostcells_p2(s_j, elem, 1);
        EnthalpyWeightedLap_bilinear_p(elem, node, s_j, hplus, hcenter, Sol, mpv, dt, t_j);
        
        omega = 0.0; 
        tmp = 0.0; 
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    omega += s_j[n] * t_j[n];
                    tmp   += t_j[n] * t_j[n];
                }
            }
        }
        
        omega /= tmp;
        
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    p2[n]    += alpha * p_j[n]  + omega * s_j[n];
                    r_j[n]    = s_j[n] - omega * t_j[n];
                }
            }
        }
        
        tmp       = 0.0;
        tmp_local = 0.0;
        for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
            for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                for(i = igx; i < icx - igx; i++) {n = m + i;
                    tmp      += r_j[n] * r_j[n];
                    tmp_local = MAX_own(tmp_local, fabs(r_j[n])); 
                }
            }
        }
        
        rho1 = rho2;
        tmp_local *= 0.5*dt/(precon_inv_scale*local_precision);
        tmp = 0.5*dt*sqrt(tmp/cell_cnt)/(precon_inv_scale*precision);
        
        cnt++;
        if(cnt % 100 == 0) printf(" iter = %d,  residual = %e\n", cnt, tmp);
    }
    
    printf(" iter = %d,  residual = %e,  local residual = %e,  gridsize = %d\n", cnt, tmp, tmp_local, nc);
    
    data->actual_iterations = cnt;
    
    return(MIN_own(tmp,tmp_local));
}

#endif /* SOLVER_1_CR2 */

