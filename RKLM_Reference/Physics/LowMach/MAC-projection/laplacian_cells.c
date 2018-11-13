/*******************************************************************************
File:   laplacian_cells.c
Author: Nicola                        Rupert     
Date:   Fri Mar 13 07:56:56 CET 1998  Feb. 2004
*******************************************************************************/

#include <assert.h> 
#include <math.h>
#include <float.h>
#include <string.h>

#include "Common.h"
#include "laplacian_cells.h"
#include "error.h"
#include "userdata.h"
#include "thermodynamic.h"
#include "EOS.h"
#include "math_own.h"
#include "ThomasAlgorithmus.h"

/* ========================================================================== */

static enum Boolean diaginv_c_is_allocated = WRONG;
static double *diaginv_c;
static double *diag_c;

/* -------------------------------------------------------------------------- */

double precon_c_diag_prepare(
                      const NodeSpaceDiscr* node,
                      const ElemSpaceDiscr* elem,
                      const double* hplus[3],
                      const double* hcenter) {
    
    /*
     the following defines diag/diaginv directly through the calculations in
     EnthalpyWeightedLap_bilinear_p()
     As a consequence, also the boundary elements have the real diagonal values,
     as opposed to what I did previously, when all the elements received the
     diagonal value from a full stencil.
     */
    
    double precon_inv_scale;
    
    if (diaginv_c_is_allocated == WRONG) {
        diag_c    = (double*)malloc(node->nc*sizeof(double));
        diaginv_c = (double*)malloc(node->nc*sizeof(double));
        diaginv_c_is_allocated = CORRECT;
    }
    
    const int ndim = elem->ndim;
    
    switch(ndim) {
        case 1: {
            const int igx       = elem->igx;
            const int icx       = elem->icx;
            const double dx     = elem->dx;
            const double oodx2  = 1.0 / (dx * dx);
            
            const double* hplusx = hplus[0];
            const double* hc     = hcenter;
                        
            for(int nc=0; nc<elem->nc; nc++) diag_c[nc] = diaginv_c[nc] = 0.0;
            
            for(int nc = igx; nc < icx - igx; nc++) {                
                
                diag_c[nc]  = - oodx2 * ( hplusx[nc+1] + hplusx[nc] );
                diag_c[nc] += hc[nc];
            }
            
            precon_inv_scale = 0.0;
            
            for (int nc = igx; nc < icx-igx; nc++) {
                diaginv_c[nc] = 1.0/diag_c[nc];
                precon_inv_scale = MAX_own(precon_inv_scale, fabs(diaginv_c[nc]));
            }
            
            break;
        }
        case 2: {
            const int igx       = elem->igx;
            const int icx       = elem->icx;
            const int ifx       = elem->ifx;
            const int igy       = elem->igy;
            const int icy       = elem->icy;
            const int ify       = elem->ify;
            const double dx     = elem->dx;
            const double dy     = elem->dy;
            const double oodx2  = 1.0 / (dx * dx);
            const double oody2  = 1.0 / (dy * dy);
                        
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
            const double* hc     = hcenter;

            const double b = P1_ALTERNATIVE_STENCIL_WEIGHT; 
            const double a = 1.0-2.0*b;
            
            for(int nc=0; nc<elem->nc; nc++) diag_c[nc] = diaginv_c[nc] = 0.0;
            
            for(int j = igy; j < icy - igy; j++) {int mc = j * icx;
                for(int i = igx; i < icx - igx; i++) {int nc  = mc + i;
                    
                    int ox  = j * ifx + i;
                    int oy  = i * ify + j;
                    int o_e = ox + 1;
                    int o_w = ox;
                    int o_n = oy + 1;
                    int o_s = oy;
                    
                    diag_c[nc]  = - a * oodx2 * ( hplusx[o_e] + hplusx[o_w] );
                    diag_c[nc] += - a * oody2 * ( hplusy[o_n] + hplusy[o_s] );
                    
                    diag_c[nc] += hc[nc];
                }
            }

            precon_inv_scale = 0.0;
            
            for (int j = igy; j < icy-igy; j++) {int mc = j*icx;
                for (int i = igx; i < icx-igx; i++) {int nc = mc + i;
                    diaginv_c[nc] = 1.0/diag_c[nc];
                    precon_inv_scale = MAX_own(precon_inv_scale, fabs(diaginv_c[nc]));
                }
            }

            break;
        }
        case 3: {
            const int igx   = elem->igx;
            const int icx   = elem->icx;
            const int ifx   = elem->ifx;
            
            const int igy   = elem->igy;
            const int icy   = elem->icy;
            const int ify   = elem->ify;
            
            const int igz   = elem->igz;
            const int icz   = elem->icz;
            const int ifz   = elem->ifz;

            const double dx = elem->dx;
            const double dy = elem->dy;
            const double dz = elem->dz;            
                        
            for(int nc=0; nc<elem->nc; nc++) diag_c[nc] = diaginv_c[nc] = 0.0;
            
            const double codsq[3] = {1.0/(64.0*dx*dx), 1.0/(64.0*dy*dy), 1.0/(64.0*dz*dz)}; 
            
            const int ics[3][3] = {{icx,icy,icz}, {icy,icz,icx}, {icz,icx,icy}};
            const int ifs[3][3] = {{ifx,ify,ifz}, {ify,ifz,ifx}, {ifz,ifx,ify}};
            
            int i, j, k, l, m, n, idim;
                        
            for(k = igz; k < icz - igz; k++) {l = k * icx*icy;                
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;                
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        
                        int ijk[3][3] = {{i,j,k}, {j,k,i}, {k,i,j}};
                        
                        for(idim = 0; idim < 3; idim++) {
                            int o_m = ijk[idim][2] * ics[idim][1] * ifs[idim][0] + ijk[idim][1] * ifs[idim][0] + ijk[idim][0];
                            int o_p = o_m + 1;
                            
                            double h_p = hplus[idim][o_p];
                            double h_m = hplus[idim][o_m];
                                                        
                            diag_c[n] -= codsq[idim] * 36.0 * (h_p + h_m);
                        }
                        diag_c[n] += hcenter[n];
                    }
                }    
            }
                        
            precon_inv_scale = 0.0;
            
            for (int k = igz; k < icz-igz; k++) {int lc = k*icx*icy;
                for (int j = igy; j < icy-igy; j++) {int mc = lc + j*icx;
                    for (int i = igx; i < icx-igx; i++) {int nc = mc + i;
                        diaginv_c[nc] = 1.0/diag_c[nc];
                        precon_inv_scale = MAX_own(precon_inv_scale, fabs(diaginv_c[nc]));
                    }
                }
            }
            
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
    
    return precon_inv_scale;

}

/* -------------------------------------------------------------------------- */

void precon_c_diag_apply(
                    double* vec_out,
                    const double* vec_in,
                    const ElemSpaceDiscr *elem) {
    for (int nc=0; nc<elem->nc; nc++) {
        vec_out[nc] = vec_in[nc]*diag_c[nc];
    }
}

/* -------------------------------------------------------------------------- */

void precon_c_diag_invert(
                     double* vec_out,
                     const double* vec_in,
                     const ElemSpaceDiscr *elem) {
    for (int nc=0; nc<elem->nc; nc++) {
        vec_out[nc] = vec_in[nc]*diaginv_c[nc];
    }
}

/* ========================================================================== */

static enum Boolean tridiago_is_allocated = WRONG;
static double *diaginv_c;
static double *diag_c;
static double *tridiago[3];
static double* upper_c;
static double* diago_c;
static double* lower_c;
static double* v_in_c;
static double* v_out_c;
static int size;

/* -------------------------------------------------------------------------- */

double precon_c_column_prepare(
                        const NodeSpaceDiscr* node,
                        const ElemSpaceDiscr* elem,
                        const double* hplus[3],
                        const double* hcenter) {
    
    /*
     the following defines diag/diaginv directly through the calculations in
     EnthalpyWeightedLap_bilinear_p()
     As a consequence, also the boundary elements have the real diagonal values,
     as opposed to what I did previously, when all the elements received the
     diagonal value from a full stencil.
     */
    
    extern User_Data ud;
    
    double precon_inv_scale, diag_scale;
        
    size = elem->icy-2*elem->igy;
    
    if (tridiago_is_allocated == WRONG) {
        for (int k=0; k<3; k++) {
            tridiago[k] = (double*)malloc(elem->nc*sizeof(double));
        }

        diag_c    = (double*)malloc(node->nc*sizeof(double));
        diaginv_c = (double*)malloc(node->nc*sizeof(double));
        upper_c   = (double*)malloc(size*sizeof(double));
        diago_c   = (double*)malloc(size*sizeof(double));
        lower_c   = (double*)malloc(size*sizeof(double));
        v_in_c    = (double*)malloc(size*sizeof(double));
        v_out_c   = (double*)malloc(size*sizeof(double));

        tridiago_is_allocated = CORRECT;
    }
    
    const int ndim = elem->ndim;
    
    switch(ndim) {
        case 1: {
            ERROR("function not available for 1D");
            break;
        }
        case 2: {
            const int igx       = elem->igx;
            const int icx       = elem->icx;
            const int ifx       = elem->ifx;
            const int igy       = elem->igy;
            const int icy       = elem->icy;
            const int ify       = elem->ify;
            const double dx     = elem->dx;
            const double dy     = elem->dy;
            const double oody2  = 1.0 / (dy * dy);
            const double oodx2  = 1.0 / (dx * dx);
                        
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
            const double* hc     = hcenter;
            
            for (int k=0; k<3; k++) {
                for(int nn=0; nn<elem->nc; nn++) tridiago[k][nn] = 0.0;
            }
            
            for(int j = igy; j < icy - igy; j++) {int mc = j * icx;
                                
                for(int i = igx; i < icx - igx; i++) {int nc  = mc + i;
                    
                    int ox = j * ifx + i;
                    int o_e   = ox + 1;
                    int o_w   = ox;

                    int oy  = i * ify + j;
                    int o_n = oy + 1;
                    int o_s = oy;
                    
                    /* !!!! this could need fixing for the  P1_ALTERNATIVE_STENCIL  option ! 
                    assert(0);
                     */
                    
                    tridiago[0][nc] = oody2 * hplusy[o_s]; 
                    tridiago[1][nc] = - oody2 * (hplusy[o_n] + hplusy[o_s])
                                      - oodx2 * (hplusx[o_e] + hplusx[o_w]) + hc[nc];
                    tridiago[2][nc] = oody2 * hplusy[o_n]; 

                    
                    diag_c[nc]  = oodx2 * ( -(hplusx[o_e] + hplusx[o_w])
                                           + 0.125 * (  hplusx[o_e] + hplusx[o_w] + hplusx[o_e] + hplusx[o_w])
                                           );
                    
                    diag_c[nc] += oody2 * ( -(hplusy[o_n] + hplusy[o_s])
                                           + 0.125 * (  hplusy[o_n] + hplusy[o_s] + hplusy[o_n] + hplusy[o_s])
                                           );
                    
                    diag_c[nc] += hc[nc];
                    
                }
            }
                        
            diag_scale = 1.0;
            
            precon_inv_scale = 0.0;
            
            for (int j = igy; j < icy-igy; j++) {int mc = j*icx;
                for (int i = igx; i < icx-igx; i++) {int nc = mc + i;
                    diag_c[nc]   *= diag_scale;
                    diaginv_c[nc] = 1.0/diag_c[nc];
                    precon_inv_scale = MAX_own(precon_inv_scale, fabs(diaginv_c[nc]));
                }
            }

            break;
        }
        case 3: {
            const int igx       = elem->igx;
            const int icx       = elem->icx;
            const int ifx       = elem->ifx;
            const int igy       = elem->igy;
            const int icy       = elem->icy;
            const int ify       = elem->ify;
            const int igz       = elem->igz;
            const int icz       = elem->icz;
            const int ifz       = elem->ifz;
            const double dx     = elem->dx;
            const double dy     = elem->dy;
            const double dz     = elem->dz;
            const double oody2  = 1.0 / (dy * dy);
            const double oodx2  = 1.0 / (dx * dx);
            const double oodz2  = 1.0 / (dz * dz);
            
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
            const double* hplusz = hplus[2];
            const double* hc     = hcenter;
            
            for (int k=0; k<3; k++) {
                for(int nn=0; nn<elem->nc; nn++) tridiago[k][nn] = 0.0;
            }
            
            for(int k = igz; k < icz - igz; k++) {
                int lc = k * icx*icy;
                
                for(int j = igy; j < icy - igy; j++) {
                    int mc = lc + j * icx;
                    
                    for(int i = igx; i < icx - igx; i++) {
                        int nc  = mc + i;
                        
                        int ox  = k*ifx*icy + j*ifx + i;
                        int o_e = ox + 1;
                        int o_w = ox;
                        
                        int oy  = i*ify*icz + k*ify + j;
                        int o_n = oy + 1;
                        int o_s = oy;
                        
                        int oz  = j*ifz*icx + i*ifz + k;
                        int o_f = oz + 1;
                        int o_b = oz;
                        
                        tridiago[0][nc] = oody2 * hplusy[o_s]; 
                        
                        tridiago[1][nc] = - oodx2 * (hplusx[o_e] + hplusx[o_w]) 
                                          - oody2 * (hplusy[o_n] + hplusy[o_s]) 
                                          - oodz2 * (hplusz[o_f] + hplusz[o_b]) 
                                          + hc[nc];
                        
                        tridiago[2][nc] = oody2 * hplusy[o_n]; 
                        
                        diag_c[nc]  = oodx2 * ( -0.5625*(hplusx[o_e] + hplusx[o_w]) );
                        
                        diag_c[nc] += oody2 * ( -0.5625*(hplusy[o_n] + hplusy[o_s]) );

                        diag_c[nc] += oodz2 * ( -0.5625*(hplusz[o_f] + hplusz[o_b]) );

                        diag_c[nc] += hc[nc];
                        
                    }
                }
            }
            
            diag_scale = 1.0;
            
            precon_inv_scale = 0.0;
            
            for (int k = igz; k < icz-igz; k++) {
                int lc = k*icx*icy;
                for (int j = igy; j < icy-igy; j++) {
                    int mc = lc + j*icx;
                    for (int i = igx; i < icx-igx; i++) {
                        int nc = mc + i;
                        diag_c[nc]   *= diag_scale;
                        diaginv_c[nc] = 1.0/diag_c[nc];
                        precon_inv_scale = MAX_own(precon_inv_scale, fabs(diaginv_c[nc]));
                    }  
                }
            }
        }
            break;
        default: 
            ERROR("ndim not in {1, 2, 3}");
            break;
    }
    
    return precon_inv_scale;
    
}

/* -------------------------------------------------------------------------- */

void precon_c_column_apply(
                    double* vec_out,
                    const double* vec_in,
                    const ElemSpaceDiscr *elem) {

    const int icxe = elem->icx;
    const int icye = elem->icy;
    const int icze = elem->icz;
    
    const int igxe = elem->igx;
    const int igye = elem->igy;
    const int igze = elem->igz;
    
    for (int k = igze; k<icze-igze; k++) {
        int lc = k*icxe*icye;
        for (int i = igxe; i<icxe-igxe; i++) {
            int mc = lc + i;
            for (int j = igye; j<icye-igye; j++) {
                int nn0 = mc + (j-1)*icxe;
                int nn1 = mc +   j  *icxe;
                int nn2 = mc + (j+1)*icxe;
                vec_out[nn1] = tridiago[0][nn1]*vec_in[nn0]+tridiago[1][nn1]*vec_in[nn1]+tridiago[2][nn1]*vec_in[nn2];
            }
        }        
    }
}

/* -------------------------------------------------------------------------- */

void precon_c_column_invert(
                     double* vec_out,
                     const double* vec_in,
                     const ElemSpaceDiscr *elem) {

    /* It will be more efficient to immediately store tridiago in the y-first ordering
     to avoid the transposition and intermediate storage in upper, diago, lower,
     but that can be done later
     */
    
    const int icxe = elem->icx;
    const int icye = elem->icy;
    const int icze = elem->icz;
    
    const int igxe = elem->igx;
    const int igye = elem->igy;
    const int igze = elem->igz;
        
    for (int k=igze; k<icze-igze; k++) {
        int lc = k*icxe*icye;
        
        for (int i=igxe; i<icxe-igxe; i++) {
            int jbot, jtop;
            int mc = lc + i;
            
            for (int j=igye; j<icye-igye; j++) {
                int j_inc = j-igye;
                int nc    = mc + j*icxe;
                lower_c[j_inc] = tridiago[0][nc]*diaginv_c[nc];
                diago_c[j_inc] = tridiago[1][nc]*diaginv_c[nc];
                upper_c[j_inc] = tridiago[2][nc]*diaginv_c[nc];
                v_in_c[j_inc]  = vec_in[nc]*diaginv_c[nc];
            }
            
            /* fix diagonal entries in the first and last row for boundary condition consistency */
            jbot = 0;
            jtop = icye-2*igye-1;
            
            diago_c[jbot] += lower_c[jbot];
            lower_c[jbot]  = 0.0;
            diago_c[jtop] += upper_c[jtop];
            upper_c[jtop]  = 0.0;
            
            Thomas_Algorithm(v_out_c, v_in_c, upper_c, diago_c, lower_c, size);
            
            for (int j=igye; j<icye-igye; j++) {
                int j_inc = j-igye;
                int nc    = mc + j*icxe;
                vec_out[nc] = v_out_c[j_inc];
            }
        }
    }
}


/* -------------------------------------------------------------------------- */

double precon_c_prepare(
                      const NodeSpaceDiscr* node,
                      const ElemSpaceDiscr* elem,
                      const double* hplus[3],
                      const double* hcenter) {

    extern User_Data ud;
    
    if (ud.gravity_strength[ud.gravity_direction] > 0 && ud.column_preconditioner == CORRECT) {
        return precon_c_column_prepare(node, elem, hplus, hcenter);
    } else {
        return precon_c_diag_prepare(node, elem, hplus, hcenter);        
    }
}

/* -------------------------------------------------------------------------- */

void precon_c_apply(
                    double* vec_out,
                    const double* vec_in,
                    const ElemSpaceDiscr *elem) 
{
#ifdef PRECON    
    /* TODO: controlled redo of changes from 2018.10.24 to 2018.11.11 
     make sure the preconditioning calls are issued for  October 24 version */

    extern User_Data ud;
    
    if (ud.gravity_strength[ud.gravity_direction] > 0 && ud.column_preconditioner == CORRECT) {
        precon_c_column_apply(vec_out, vec_in, elem);
    } else {
        precon_c_diag_apply(vec_out, vec_in, elem);        
    }
#else
    return;
#endif
    
}

/* -------------------------------------------------------------------------- */

void precon_c_invert(
                     double* vec_out,
                     const double* vec_in,
                     const ElemSpaceDiscr *elem) 
{
#ifdef PRECON
    /* TODO: controlled redo of changes from 2018.10.24 to 2018.11.11 
     make sure the preconditioning calls are issued for  October 24 version */

    extern User_Data ud;
    
    if (ud.gravity_strength[ud.gravity_direction] > 0 && ud.column_preconditioner == CORRECT) {
        precon_c_column_invert(vec_out, vec_in, elem);
    } else {
        precon_c_diag_invert(vec_out, vec_in, elem);        
    }
#else
    return;
#endif
}


/* ========================================================================== */
#define OUTPUT_LAP_CELLS 0
#if OUTPUT_LAP_CELLS
#include "io.h"
static int lap_output_count = 0;
#endif

void EnthalpyWeightedLap_bilinear_p(
                                    const ElemSpaceDiscr* elem, 
                                    const NodeSpaceDiscr* node,
                                    const double* p,
                                    const double* hplus[3],
                                    const double* hcenter,
                                    const MPV* mpv, 
                                    double* lap) {
    
    const int ndim = elem->ndim;
    
    switch(ndim) {
        case 1: {
            const int igx = elem->igx;
            const int icx = elem->icx;
            
            const double dx    = elem->dx;            
            const double oodx2 = 1.0 / (dx * dx);
            
            const double* hplusx = hplus[0];
            const double* hc     = hcenter;
                                                
            for(int nc = igx; nc < icx - igx; nc++) {
                lap[nc]  = oodx2 * ( hplusx[nc+1] * (p[nc+1] - p[nc]) - hplusx[nc] * (p[nc] - p[nc-1]) );
                lap[nc] += hc[nc] * p[nc];
            }
            precon_c_invert(lap, lap, elem);
            break;
        }
        case 2: {
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            const double dx = elem->dx;
            const double dy = elem->dy;
                    
#ifdef NEW_LAP
            const double oodx    = 1.0 / dx;
            const double oody    = 1.0 / dy;
#else
            const double oodx2   = 1.0 / (dx * dx);
            const double oody2   = 1.0 / (dy * dy);
#endif       
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
            const double* hc     = hcenter;
            
            int n_c, n_e, n_ne, n_n, n_nw, n_w, n_sw, n_s, n_se;
            int o_e, o_n, o_w, o_s, ox, oy;
            
            int i, j, m, n;
            
            const double b = P1_ALTERNATIVE_STENCIL_WEIGHT;  
            const double a = 1.0-2.0*b;

            
            for(j = igy; j < icy - igy; j++) {m = j * icx;				
                for(i = igx; i < icx - igx; i++) {n  = m + i;
                    
                    ox = j * ifx + i;
                    oy = i * ify + j;
                    
                    n_c     = n;
                    n_e     = n + 1      ;
                    n_ne    = n + 1 + icx;
                    n_n     = n     + icx;
                    n_nw    = n - 1 + icx;
                    n_w     = n - 1      ;
                    n_sw    = n - 1 - icx;
                    n_s     = n     - icx;
                    n_se    = n + 1 - icx;
                    
                    o_e   = ox + 1;
                    o_w   = ox;
                    
                    o_n   = oy + 1;
                    o_s   = oy;
                    
#ifdef NEW_LAP
                    double dpdx_e  = oodx * (p[n_e ] - p[n_c ]);
                    double dpdx_ne = oodx * (p[n_ne] - p[n_n ]);
                    double dpdx_se = oodx * (p[n_se] - p[n_s ]);
                    double dpdx_w  = oodx * (p[n_c ] - p[n_w ]);
                    double dpdx_nw = oodx * (p[n_n ] - p[n_nw]);
                    double dpdx_sw = oodx * (p[n_s ] - p[n_sw]);
                    
                    lap[n]  = oodx * (a * ( hplusx[o_e] * dpdx_e            - hplusx[o_w] * dpdx_w ) 
                                    + b * ( hplusx[o_e] * (dpdx_ne+dpdx_se) - hplusx[o_w] * (dpdx_nw+dpdx_sw)));
                    
                    double dpdy_n  = oody * (p[n_n ] - p[n_c ]);
                    double dpdy_ne = oody * (p[n_ne] - p[n_e ]);
                    double dpdy_nw = oody * (p[n_nw] - p[n_w ]);
                    double dpdy_s  = oody * (p[n_c ] - p[n_s ]);
                    double dpdy_se = oody * (p[n_e ] - p[n_se]);
                    double dpdy_sw = oody * (p[n_w ] - p[n_sw]);
                    
                    lap[n] += oody * (a * ( hplusy[o_n] * dpdy_n            - hplusy[o_s] * dpdy_s )
                                    + b * ( hplusy[o_n] * (dpdy_ne+dpdy_nw) - hplusy[o_s] * (dpdy_se+dpdy_sw)));
#else               
                    lap[n]  = oodx2 * (  a * ( hplusx[o_e] * (p[n_e ] - p[n_c ]) - hplusx[o_w] * (p[n_c ] - p[n_w ]) ) 
                                       + b * ( hplusx[o_e] * ( ( p[n_ne] - p[n_n ] ) + ( p[n_se] - p[n_s ] ) ) 
                                             - hplusx[o_w] * ( ( p[n_n ] - p[n_nw] ) + ( p[n_s ] - p[n_sw] ) )  
                                             )
                                       );
                    
                    lap[n] += oody2 * (  a * ( hplusy[o_n] * (p[n_n ] - p[n_c ]) - hplusy[o_s] * (p[n_c ] - p[n_s ]) )
                                       + b * ( hplusy[o_n] * ( ( p[n_ne] - p[n_e ] ) + ( p[n_nw] - p[n_w ] ) ) 
                                             - hplusy[o_s] * ( ( p[n_e ] - p[n_se] ) + ( p[n_w ] - p[n_sw] ) ) 
                                             )
                                       );
#endif
                    
                    lap[n] += hc[n] * p[n_c];
                    
                }
            }    
                        
            precon_c_invert(lap, lap, elem);

            break;
        }
        case 3: {
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            const int igz = elem->igz;
            const int icz = elem->icz;
            const int ifz = elem->ifz;
            
            const double dx = elem->dx;
            const double dy = elem->dy;
            const double dz = elem->dz;
            
#ifdef CORIOLIS_EXPLICIT
            /* const double coriolis  = 0.0; */
#else
            assert(0); /* Modification of Laplacian for implicit Coriolis not implemented yet */
            /* const double coriolis  = ud.coriolis_strength[0]; */
#endif

            /* for MG-scaling with elem->scale_factor; see old version of routine */
            
            const double codsq[3] = {1.0/(64.0*dx*dx), 1.0/(64.0*dy*dy), 1.0/(64.0*dz*dz)}; 
            
            const int dis[3][3] = {{1, icx, icx*icy}, {icx, icx*icy, 1}, {icx*icy, 1, icx}};
            const int ics[3][3] = {{icx,icy,icz}, {icy,icz,icx}, {icz,icx,icy}};
            const int ifs[3][3] = {{ifx,ify,ifz}, {ify,ifz,ifx}, {ifz,ifx,ify}};
            
            int i, j, k, l, m, n, idim;
            
            memset(lap, 0.0, elem->nc*sizeof(double));
            
            for(k = igz; k < icz - igz; k++) {l = k * icx*icy;				
                for(j = igy; j < icy - igy; j++) {m = l + j * icx;				
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        
                        int ijk[3][3] = {{i,j,k}, {j,k,i}, {k,i,j}};
                        
                        for(idim = 0; idim < 3; idim++) {
                            int n_m     = n - dis[idim][0];
                            int n_c     = n;
                            int n_p     = n + dis[idim][0];
                            
                            int o_m = ijk[idim][2] * ics[idim][1] * ifs[idim][0] + ijk[idim][1] * ifs[idim][0] + ijk[idim][0];
                            int o_p = o_m + 1;
                            
                            double h_p = hplus[idim][o_p];
                            double h_m = hplus[idim][o_m];
                            
                            int n_ms    = n - dis[idim][0] - dis[idim][1];
                            int n_cs    = n                - dis[idim][1];
                            int n_ps    = n + dis[idim][0] - dis[idim][1];
                            
                            int n_me    = n - dis[idim][0] + dis[idim][2];
                            int n_ce    = n                + dis[idim][2];
                            int n_pe    = n + dis[idim][0] + dis[idim][2];
                            
                            int n_mn    = n - dis[idim][0] + dis[idim][1];
                            int n_cn    = n                + dis[idim][1];
                            int n_pn    = n + dis[idim][0] + dis[idim][1];
                            
                            int n_mw    = n - dis[idim][0] - dis[idim][2];
                            int n_cw    = n                - dis[idim][2];
                            int n_pw    = n + dis[idim][0] - dis[idim][2];
                            
                            int n_msw   = n - dis[idim][0] - dis[idim][1] - dis[idim][2];
                            int n_csw   = n                - dis[idim][1] - dis[idim][2];
                            int n_psw   = n + dis[idim][0] - dis[idim][1] - dis[idim][2];
                            
                            int n_mse   = n - dis[idim][0] - dis[idim][1] + dis[idim][2];
                            int n_cse   = n                - dis[idim][1] + dis[idim][2];
                            int n_pse   = n + dis[idim][0] - dis[idim][1] + dis[idim][2];
                            
                            int n_mne   = n - dis[idim][0] + dis[idim][1] + dis[idim][2];
                            int n_cne   = n                + dis[idim][1] + dis[idim][2];
                            int n_pne   = n + dis[idim][0] + dis[idim][1] + dis[idim][2];
                            
                            int n_mnw   = n - dis[idim][0] + dis[idim][1] - dis[idim][2];
                            int n_cnw   = n                + dis[idim][1] - dis[idim][2];
                            int n_pnw   = n + dis[idim][0] + dis[idim][1] - dis[idim][2];
                            
                            lap[n] += codsq[idim] *                                            \
                            (  36.0 * (  h_p * (p[n_p] - p[n_c]) - h_m * (p[n_c] - p[n_m])          \
                                       )                                                            \
                             +  6.0 * (  h_p * (p[n_ps] - p[n_cs]) - h_m * (p[n_cs] - p[n_ms])      \
                                       + h_p * (p[n_pe] - p[n_ce]) - h_m * (p[n_ce] - p[n_me])      \
                                       + h_p * (p[n_pn] - p[n_cn]) - h_m * (p[n_cn] - p[n_mn])      \
                                       + h_p * (p[n_pw] - p[n_cw]) - h_m * (p[n_cw] - p[n_mw])      \
                                       )                                                            \
                             +  1.0 * (  h_p * (p[n_psw] - p[n_csw]) - h_m * (p[n_csw] - p[n_msw])  \
                                       + h_p * (p[n_pse] - p[n_cse]) - h_m * (p[n_cse] - p[n_mse])  \
                                       + h_p * (p[n_pne] - p[n_cne]) - h_m * (p[n_cne] - p[n_mne])  \
                                       + h_p * (p[n_pnw] - p[n_cnw]) - h_m * (p[n_cnw] - p[n_mnw])  \
                                       )                                                            \
                             );
                        }
                        lap[n] += hcenter[n] * p[n];
                    }
                }    
            }    
            
            precon_c_invert(lap, lap, elem);

            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
    
#if OUTPUT_LAP_CELLS
    extern User_Data ud;
    FILE *plapfile = NULL;
    char fn[120], fieldname[90];
    if (lap_output_count < 10) {
        sprintf(fn, "%s/lap_cells/lap_cells_00%d.hdf", ud.file_name, lap_output_count);
    } else if(lap_output_count < 100) {
        sprintf(fn, "%s/lap_cells/lap_cells_0%d.hdf", ud.file_name, lap_output_count);
    } else {
        sprintf(fn, "%s/lap_cells/lap_cells_%d.hdf", ud.file_name, lap_output_count);
    }
    sprintf(fieldname, "lap_cells");    
    WriteHDF(plapfile, elem->icx, elem->icy, elem->icz, elem->ndim, lap, fn, fieldname);
    
    lap_output_count++;
#endif

}
