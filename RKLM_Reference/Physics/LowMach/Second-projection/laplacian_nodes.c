/*******************************************************************************
 File:   laplacian_nodes.c
 Author: Nicola                        Rupert
 Date:   Fri Mar 13 07:56:56 CET 1998  Feb. 2004
 *******************************************************************************/
#include <float.h>
#include <string.h>
#include <math.h>
#include "laplacian_nodes.h"
#include "userdata.h"
#include "error.h"
#include "Common.h"
#include "enumerator.h"
#include "ThomasAlgorithmus.h"
#include "mpv.h"
#include "math_own.h"


static enum Boolean diaginv_is_allocated = WRONG;
static double *diaginv;
static double *diag;

static enum Boolean tridiago_is_allocated = WRONG;
static double *tridiago[3];

/* ========================================================================== */

double precon_diag_prepare(
                    const NodeSpaceDiscr* node,
                    const ElemSpaceDiscr* elem,
                    const double* hplus[3],
                    const double* wcenter,
                    const int x_periodic,
                    const int y_periodic,
                    const int z_periodic) 
{
    
    /*
     the following defines diag/diaginv directly through the calculations in
     EnthalpyWeightedLap_Node_bilinear_p_scatter()
     As a consequence, also the boundary elements have the real diagonal values,
     as opposed to what I did previously, when all the elements received the
     diagonal value from a full stencil.
     */
    
    extern User_Data ud;
    extern MPV* mpv;
    
    double precon_inv_scale;
    
    if (diaginv_is_allocated == WRONG) {
        diag    = (double*)malloc(node->nc*sizeof(double));
        diaginv = (double*)malloc(node->nc*sizeof(double));
        diaginv_is_allocated = CORRECT;
    }
    
    const int ndim = node->ndim;
    
    switch(ndim) {
        case 1: {
            ERROR("function not available");
            break;
        }
        case 2: {
            const int igxn = node->igx;
            const int icxn = node->icx;
            const int igyn = node->igy;
            const int icyn = node->icy;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            
            const double dx = node->dx;
            const double dy = node->dy;
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            const double* hcenter  = wcenter;
                        
            const double oodx2 = 0.5 / (dx * dx);
            const double oody2 = 0.5 / (dy * dy);
            const double nine_pt = (0.25 * (1.0 + P2_DIAGONAL_FIVE_POINT)) * P2_FULL_STENCIL;
            
            double flux_x_lower, flux_x_upper, flux_y_left, flux_y_right, hc;
            double dsq_p_dxdy;
            
            double pnn      = 1.0;
            double pnn1     = 1.0;
            double pnnicxn  = 1.0;
            double pnn1icxn = 1.0;
            
            int i, j, me, mn, ne, nn, nn1, nnicxn, nn1icxn;
            
            for(nn=0; nn<node->nc; nn++) diag[nn] = diaginv[nn] = 0.0;
            
            for(j = igye; j < icye - igye; j++) {
                me   = j * icxe;
                mn   = j * icxn;
                
                for(i = igxe; i < icxe - igxe; i++) {
                    ne       = me + i;
                    
                    nn       = mn + i;
                    nn1      = nn + 1;
                    nnicxn   = nn + icxn;
                    nn1icxn  = nn + 1 + icxn;
                    
                    /* location nn */
                    dsq_p_dxdy    = pnn;
                    
                    flux_x_lower  = hplusx[ne] * oodx2 * ( - pnn + nine_pt * dsq_p_dxdy);
                    flux_x_upper  = hplusx[ne] * oodx2 * (       - nine_pt * dsq_p_dxdy);
                    
                    flux_y_left   = hplusy[ne] * oody2 * ( - pnn + nine_pt * dsq_p_dxdy);
                    flux_y_right  = hplusy[ne] * oody2 * (       - nine_pt * dsq_p_dxdy);
                    
                    hc            = 0.25 * hcenter[ne];
                    
                    diag[nn]      += (  flux_x_lower + flux_y_left ) + hc*pnn;
                    
                    
                    
                    /* location nn1 */
                    dsq_p_dxdy    = - pnn1;
                    
                    flux_x_lower  = hplusx[ne] * oodx2 * ( pnn1 + nine_pt * dsq_p_dxdy);
                    flux_x_upper  = hplusx[ne] * oodx2 * (      - nine_pt * dsq_p_dxdy);
                    
                    flux_y_left   = hplusy[ne] * oody2 * (      + nine_pt * dsq_p_dxdy);
                    flux_y_right  = hplusy[ne] * oody2 * (-pnn1 - nine_pt * dsq_p_dxdy);
                    
                    hc            = 0.25 * hcenter[ne];
                    
                    diag[nn1]     += (- flux_x_lower + flux_y_right) + hc*pnn1;
                    
                    
                    
                    /* location nnicxn */
                    dsq_p_dxdy    = - pnnicxn;
                    
                    flux_x_lower  = hplusx[ne] * oodx2 * (         + nine_pt * dsq_p_dxdy);
                    flux_x_upper  = hplusx[ne] * oodx2 * (-pnnicxn - nine_pt * dsq_p_dxdy);
                    
                    flux_y_left   = hplusy[ne] * oody2 * ( pnnicxn + nine_pt * dsq_p_dxdy);
                    flux_y_right  = hplusy[ne] * oody2 * (         - nine_pt * dsq_p_dxdy);
                    
                    hc            = 0.25 * hcenter[ne];
                    
                    diag[nnicxn]  += (  flux_x_upper - flux_y_left ) + hc*pnnicxn ;
                    
                    
                    
                    /* location nn1icxn */
                    dsq_p_dxdy    = pnn1icxn;
                    
                    flux_x_lower  = hplusx[ne] * oodx2 * (          + nine_pt * dsq_p_dxdy);
                    flux_x_upper  = hplusx[ne] * oodx2 * ( pnn1icxn - nine_pt * dsq_p_dxdy);
                    
                    flux_y_left   = hplusy[ne] * oody2 * (          + nine_pt * dsq_p_dxdy);
                    flux_y_right  = hplusy[ne] * oody2 * ( pnn1icxn - nine_pt * dsq_p_dxdy);
                    
                    hc            = 0.25 * hcenter[ne];
                    
                    diag[nn1icxn] += (- flux_x_upper - flux_y_right) + hc*pnn1icxn;
                }
            }
            
            /* periodicity */
            if (x_periodic) {
                for(j=igyn; j<icyn-igyn; j++) {
                    int nleft  = j * icxn + igxn;
                    int nright = j * icxn + icxn - igxn - 1;
                    
                    diag[nleft]  += diag[nright];
                    diag[nright]  = 9999.0;
                }
            }
            
            if (y_periodic) {
                for(i=igxn; i<icxn-igxn; i++) {
                    int nbottom  = i + igyn * icxn;
                    int ntop     = i + (icyn - igyn - 1) * icxn;
                    
                    diag[nbottom]  += diag[ntop];
                    diag[ntop]      = 7777.0;
                }
            }
             
            
            /* prepare for the inversion of the preconditioner */
            precon_inv_scale = 0.0;
            
            for (j = igyn; j < icyn-igyn; j++) {mn = j*icxn;
                for (i = igxn; i < icxn-igxn; i++) {nn = mn + i;
                    diaginv[nn] = 1.0/diag[nn];
                    precon_inv_scale = MAX_own(precon_inv_scale, fabs(diaginv[nn]));
                }
            }
            
            break;
        }
        case 3: {
            ERROR("function not available");
            break;
        }
            
        default: ERROR("ndim not in {1, 2, 3}");
    }  
    return precon_inv_scale;
}

/* ========================================================================== */

void precon_diag_apply(
                       double* vec_out,
                       const double* vec_in,
                       const NodeSpaceDiscr *node,
                       const int x_periodic,
                       const int y_periodic,
                       const int z_periodic) {
    
    const int icx = node->icx;
    const int icy = node->icy;
    const int icz = node->icz;
    
    const int igx = node->igx;
    const int igy = node->igy;
    const int igz = node->igz;
    
    for (int k=igz; k<MAX_own(1,icz-igz-z_periodic); k++) {
        int ln = k*icx*icy;
        for (int j=igy; j<MAX_own(1, icy-igy-y_periodic); j++) {
            int mn = ln + j*icx;
            for (int i=igx; i<MAX_own(1, icx-igx-x_periodic); i++) {
                int nn = mn + i;
                vec_out[nn] = vec_in[nn]*diag[nn];
            }
        }
    }
}

/* ========================================================================== */

void precon_diag_invert(
                        double* vec_out,
                        const double* vec_in,
                        const NodeSpaceDiscr *node,
                        const int x_periodic,
                        const int y_periodic,
                        const int z_periodic) {
    
    const int icx = node->icx;
    const int icy = node->icy;
    const int icz = node->icz;
    
    const int igx = node->igx;
    const int igy = node->igy;
    const int igz = node->igz;
    
    for (int k=igz; k<MAX_own(1,icz-igz-z_periodic); k++) {
        int ln = k*icx*icy;
        for (int j=igy; j<MAX_own(1, icy-igy-y_periodic); j++) {
            int mn = ln + j*icx;
            for (int i=igx; i<MAX_own(1, icx-igx-x_periodic); i++) {
                int nn = mn + i;
                vec_out[nn] = vec_in[nn]*diaginv[nn];
            }
        }
    }
}

/* ========================================================================== */

double precon_column_prepare(
                    const NodeSpaceDiscr* node,
                    const ElemSpaceDiscr* elem,
                    const double* hplus[3],
                    const double* wcenter,
                    const int x_periodic,
                    const int y_periodic,
                    const int z_periodic) {
    
    /*
     the following defines diag/diaginv directly through the calculations in
     EnthalpyWeightedLap_Node_bilinear_p_scatter()
     As a consequence, also the boundary elements have the real diagonal values,
     as opposed to what I did previously, when all the elements received the
     diagonal value from a full stencil.
     */
    
    double *tridiaux[3];
    
    double precon_inv_scale;
    
    if (tridiago_is_allocated == WRONG) {
        for (int k=0; k<3; k++) {
            tridiago[k] = (double*)malloc(node->nc*sizeof(double));
            tridiaux[k] = tridiago[k];
        }
        tridiago_is_allocated = CORRECT;
    }
    
    const int ndim = node->ndim;
    
    switch(ndim) {
        case 1: {
            ERROR("function not available");
            break;
        }
        case 2: {
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;

            const int icxn = node->icx;
            const int icyn = node->icy;
            const int igxn = node->igx;
            const int igyn = node->igy;

            const double dx = node->dx;
            const double dy = node->dy;
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            const double* hcenter  = wcenter;
                        
            const double oodx2 = 0.5 / (dx * dx);
            const double oody2 = 0.5 / (dy * dy);
            
            double flux_x_lower, flux_x_upper, flux_y_left, flux_y_right, hc;
            
            int i, j, me, mn, ne, nn, nn1, nnicxn, nn1icxn;
                        
            for (int k=0; k<3; k++) {
                for(nn=0; nn<node->nc; nn++) tridiago[k][nn] = 0.0;
            }
            
            for(j = igye; j < icye - igye; j++) {
                me   = j * icxe;
                mn   = j * icxn;
                
                for(i = igxe; i < icxe - igxe; i++) {
                    ne       = me + i;
                    
                    nn       = mn + i;
                    nn1      = nn + 1;
                    nnicxn   = nn + icxn;
                    nn1icxn  = nn + 1 + icxn;
                    
                    /*
                     flux_y_left   = hplusy[ne] * oody2 * ( (pnnicxn  - pnn    ) );
                     flux_y_right  = hplusy[ne] * oody2 * ( (pnn1icxn - pnn1   ) );
                     */
                    flux_y_left   = hplusy[ne] * oody2;
                    flux_y_right  = hplusy[ne] * oody2;
                    
                    flux_x_lower  = hplusx[ne] * oodx2;
                    flux_x_upper  = hplusx[ne] * oodx2;
                    
                    
                    /* summand 1.0 takes care of singularity of the tridiagoonal matrix */
                    hc            = 0.25 * hcenter[ne];
                                        
                    
                    /* eventually I should transpose the tridiago[][] field for better memory efficiency */
                    /* eventually I should let tridiago[] run from  -1 to +1 */
                    
                    tridiago[1][nn]      += - flux_x_lower - flux_y_left + hc;
                    tridiago[2][nn]      +=   flux_y_left;
                    
                    tridiago[1][nn1]     += - flux_x_lower - flux_y_right + hc;
                    tridiago[2][nn1]     +=   flux_y_right;
                    
                    tridiago[0][nnicxn]  +=   flux_y_left;
                    tridiago[1][nnicxn]  += - flux_x_upper - flux_y_left + hc;
                    
                    tridiago[0][nn1icxn] +=   flux_y_right;
                    tridiago[1][nn1icxn] += - flux_x_upper - flux_y_right + hc;                    
                }
            }
            
            /* normalization */
            precon_inv_scale = 0.0;

            for(j = igyn; j < icyn - igyn; j++) {
                mn   = j * icxn;
                
                for(i = igxn; i < icxn - igxn; i++) {
                    nn       = mn + i;

                    precon_inv_scale = MAX_own(precon_inv_scale, fabs(1.0/tridiago[1][nn]));                
                }
            }            
        }
            break;
        case 3: {
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            const int igze = elem->igz;
            const int icze = elem->icz;
            
            const int icxn = node->icx;
            const int igxn = node->igx;
            const int icyn = node->icy;
            const int igyn = node->igy;
            const int iczn = node->icz;
            const int igzn = node->igz;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            const double* hplusz   = hplus[2];
            const double* hcenter  = wcenter;
            
            const double oodx2 = 0.5 / (dx * dx);
            const double oody2 = 0.5 / (dy * dy);
            const double oodz2 = 0.5 / (dz * dz);
            
            double flux_x, flux_y, flux_z;
            double hc;
            
            int i, j, k, le, ln, me, mn, ne, nn; 
            int nn000, nn001, nn010, nn011;
            int nn100, nn101, nn110, nn111;
            
            for (int k=0; k<3; k++) {
                for(nn=0; nn<node->nc; nn++) tridiago[k][nn] = 0.0;
            }
            
            for (k = igze; k < icze - igze; k++) {
                le = k * icxe*icye;
                ln = k * icxn*icyn;
                
                for(j = igye; j < icye - igye; j++) {
                    me   = le + j * icxe;
                    mn   = ln + j * icxn;
                    
                    for(i = igxe; i < icxe - igxe; i++) {
                        ne    = me + i;
                        nn    = mn + i; 
                        
                        nn000 = nn;
                        nn001 = nn + 1;
                        nn010 = nn + icxn;
                        nn011 = nn + 1 + icxn;
                        nn100 = nn + icxn*icyn;
                        nn101 = nn + icxn*icyn + 1;
                        nn110 = nn + icxn*icyn + icxn;
                        nn111 = nn + icxn*icyn + icxn + 1;
                        
                        flux_x  = hplusx[ne] * oodx2;
                        flux_y  = hplusy[ne] * oody2;
                        flux_z  = hplusz[ne] * oodz2;
                                    
                        
                        /* summand 1.0 takes care of singularity of the tridiagoonal matrix */
                        hc            = 0.25 * hcenter[ne];
                        
                        /* eventually I should transpose the tridiago[][] field for better memory efficiency */
                        /* eventually I should let tridiago[] run from  -1 to +1 */
                        
                        tridiago[1][nn000] += - flux_x - flux_y - flux_z + hc;
                        tridiago[2][nn000] +=   flux_y;
                        
                        tridiago[1][nn001] += - flux_x - flux_y - flux_z + hc;
                        tridiago[2][nn001] +=   flux_y;
                        
                        tridiago[0][nn010] +=   flux_y;
                        tridiago[1][nn010] += - flux_x - flux_y - flux_z + hc;
                        
                        tridiago[0][nn011] +=   flux_y;
                        tridiago[1][nn011] += - flux_x - flux_y - flux_z + hc;                    

                        tridiago[1][nn100] += - flux_x - flux_y - flux_z + hc;
                        tridiago[2][nn100] +=   flux_y;
                        
                        tridiago[1][nn101] += - flux_x - flux_y - flux_z + hc;
                        tridiago[2][nn101] +=   flux_y;
                        
                        tridiago[0][nn110] +=   flux_y;
                        tridiago[1][nn110] += - flux_x - flux_y - flux_z + hc;
                        
                        tridiago[0][nn111] +=   flux_y;
                        tridiago[1][nn111] += - flux_x - flux_y - flux_z + hc;                    
                    }
                }
            }
            
            /* normalization */
            precon_inv_scale = 0.0;
            
            for(k = igzn; k < iczn - igzn; k++) {
                ln   = k*icxn*icyn;
                for(j = igyn; j < icyn - igyn; j++) {
                    mn   = ln + j*icxn;
                    for(i = igxn; i < icxn - igxn; i++) {
                        nn = mn + i;
                        precon_inv_scale = MAX_own(precon_inv_scale, fabs(1.0/tridiago[1][nn]));                
                    }
                }
            }
        }
            break;
            
        default: 
            ERROR("ndim not in {1, 2, 3}");
    }

    return precon_inv_scale;
}

/* ========================================================================== */

void precon_column_apply(
                  double* vec_out,
                  const double* vec_in,
                         const NodeSpaceDiscr *node,
                         const int x_periodic,
                         const int y_periodic,
                         const int z_periodic) {
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    const int iczn = node->icz;
    
    const int igxn = node->igx;
    const int igyn = node->igy;
    const int igzn = node->igz;
    
    for (int k = igzn; k<iczn-igzn; k++) {
        int ln = k*icxn*icyn;
        for (int i = igxn; i<icxn-igxn; i++) {
            int mn = ln + i;
            for (int j = igyn; j<icyn-igyn; j++) {
                int nn0 = mn + (j-1)*icxn;
                int nn1 = mn +   j  *icxn;
                int nn2 = mn + (j+1)*icxn;
                vec_out[nn1] = tridiago[0][nn1]*vec_in[nn0]+tridiago[1][nn1]*vec_in[nn1]+tridiago[2][nn1]*vec_in[nn2];
            }
        }
    }
}

/* ========================================================================== */

void precon_column_invert(
                   double* vec_out,
                   const double* vec_in,
                          const NodeSpaceDiscr *node,
                          const int x_periodic,
                          const int y_periodic,
                          const int z_periodic) {
    
    /* It will be more efficient to immediately store tridiago in the y-first ordering
     to avoid the transposition and intermediate storage in upper, diago, lower,
     but that can be done later
     */
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    const int iczn = node->icz;
    
    const int igxn = node->igx;
    const int igyn = node->igy;
    const int igzn = node->igz;
    
    const int size = icyn-2*igyn;
    
    double* upper  = (double*)malloc(size*sizeof(double));
    double* diago  = (double*)malloc(size*sizeof(double));
    double* lower  = (double*)malloc(size*sizeof(double));
    double* v_in   = (double*)malloc(size*sizeof(double));
    double* v_out  = (double*)malloc(size*sizeof(double));
        
    for (int k=igzn; k<iczn-igzn; k++) {
        int ln = k*icxn*icyn;
        for (int i=igxn; i<icxn-igxn; i++) {
            int mn = ln + i;
            for (int j=igyn; j<icyn-igyn; j++) {
                int j_inn = j-igyn;
                int nn    = mn + j*icxn;
                lower[j_inn] = tridiago[0][nn]/tridiago[1][nn];
                diago[j_inn] = 1.0;
                upper[j_inn] = tridiago[2][nn]/tridiago[1][nn];
                v_in[j_inn]  = vec_in[nn]/tridiago[1][nn];
            }
            Thomas_Algorithm(v_out, v_in, upper, diago, lower, size);
            for (int j=igyn; j<icyn-igyn; j++) {
                int j_inn = j-igyn;
                int nn    = mn + j*icxn;
                vec_out[nn] = v_out[j_inn];
            }
        }
    }

    free(upper);
    free(diago);
    free(lower);
    free(v_in );
    free(v_out);
    
}

/* ========================================================================== */

double precon_prepare(
                      const NodeSpaceDiscr* node,
                      const ElemSpaceDiscr* elem,
                      const double* hplus[3],
                      const double* wcenter,
                      const int x_periodic,
                      const int y_periodic,
                      const int z_periodic) 
{
    extern User_Data ud;
    
    if (ud.gravity_strength[ud.gravity_direction] > 0) {
        return precon_column_prepare(node, elem, hplus, wcenter, x_periodic, y_periodic, z_periodic);
    } else {
        return precon_diag_prepare(node, elem, hplus, wcenter, x_periodic, y_periodic, z_periodic);        
    }
}

/* ========================================================================== */

void precon_apply(
                  double* vec_out,
                  const double* vec_in,
                  const NodeSpaceDiscr *node,
                  const int x_periodic,
                  const int y_periodic,
                  const int z_periodic) 
{
    extern User_Data ud;
    
    if (ud.gravity_strength[ud.gravity_direction] > 0) {
        precon_column_apply(vec_out, vec_in, node, x_periodic, y_periodic, z_periodic);
    } else {
        precon_diag_apply(vec_out, vec_in, node, x_periodic, y_periodic, z_periodic);        
    }    
}

/* ========================================================================== */

void precon_invert(
                  double* vec_out,
                  const double* vec_in,
                   const NodeSpaceDiscr *node,
                   const int x_periodic,
                   const int y_periodic,
                   const int z_periodic) 
{
    extern User_Data ud;
    
    if (ud.gravity_strength[ud.gravity_direction] > 0) {
        precon_column_invert(vec_out, vec_in, node, x_periodic, y_periodic, z_periodic);
    } else {
        precon_diag_invert(vec_out, vec_in, node, x_periodic, y_periodic, z_periodic);        
    }    
}

/* ========================================================================== */

void EnthalpyWeightedLap_Node_bilinear_p_scatter(
                                                 const NodeSpaceDiscr* node,
                                                 const ElemSpaceDiscr* elem,
                                                 const double* p,
                                                 const double* hplus[3],
                                                 const double* wcenter,
                                                 const int x_periodic,
                                                 const int y_periodic,
                                                 const int z_periodic,
                                                 double* lap) {
    const int ndim = node->ndim;
    
    switch(ndim) {
        case 1: {
            ERROR("function not available");
            break;
        }
        case 2: {
            const int igxn = node->igx;
            const int icxn = node->icx;
            const int igyn = node->igy;
            const int icyn = node->icy;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            
            const double dx = node->dx;
            const double dy = node->dy;
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            const double* hcenter  = wcenter;
                        
            const double oodx2 = 0.5 / (dx * dx);
            const double oody2 = 0.5 / (dy * dy);
            const double nine_pt = (0.25 * (1.0 + P2_DIAGONAL_FIVE_POINT)) * P2_FULL_STENCIL;
            
            double flux_x_lower, flux_x_upper, flux_y_left, flux_y_right, hc;
            double dsq_p_dxdy;
            
            int i, j, me, mn, ne, nn, nn1, nnicxn, nn1icxn;
            
            for(nn=0; nn<node->nc; nn++) lap[nn] = 0.0;
            
            for(j = igye; j < icye - igye; j++) {
                me   = j * icxe;
                mn   = j * icxn;
                
                for(i = igxe; i < icxe - igxe; i++) {
                    ne       = me + i;
                    
                    nn       = mn + i;
                    nn1      = nn + 1;
                    nnicxn   = nn + icxn;
                    nn1icxn  = nn + 1 + icxn;
                    
                    dsq_p_dxdy    = p[nn1icxn] - p[nnicxn] - p[nn1] + p[nn];
                    
                    flux_x_lower  = hplusx[ne] * oodx2 * ( (p[nn1]     - p[nn]    ) + nine_pt * dsq_p_dxdy);
                    flux_x_upper  = hplusx[ne] * oodx2 * ( (p[nn1icxn] - p[nnicxn]) - nine_pt * dsq_p_dxdy);
                    
                    flux_y_left   = hplusy[ne] * oody2 * ( (p[nnicxn]  - p[nn]    ) + nine_pt * dsq_p_dxdy);
                    flux_y_right  = hplusy[ne] * oody2 * ( (p[nn1icxn] - p[nn1]   ) - nine_pt * dsq_p_dxdy);
                    
                    hc            = 0.25 * hcenter[ne];
                    
                    lap[nn]      += (  flux_x_lower + flux_y_left ) + hc*p[nn]     ;
                    lap[nn1]     += (- flux_x_lower + flux_y_right) + hc*p[nn1]    ;
                    lap[nnicxn]  += (  flux_x_upper - flux_y_left ) + hc*p[nnicxn] ;
                    lap[nn1icxn] += (- flux_x_upper - flux_y_right) + hc*p[nn1icxn];
                }
            }
                        
            if (x_periodic) {
                for(j=igyn; j<icyn-igyn; j++) {
                    int nleft  = j * icxn + igxn;
                    int nright = j * icxn + icxn - igxn - 1;
                    
                    lap[nleft]  += lap[nright];
                    lap[nright]  = 0.0;
                }
            }
            
            if (y_periodic) {
                for(i=igxn; i<icxn-igxn; i++) {
                    int nbottom  = i + igyn * icxn;
                    int ntop     = i + (icyn - igyn - 1) * icxn;
                    
                    lap[nbottom]  += lap[ntop];
                    lap[ntop]      = 0.0;
                }
            }
            
            precon_invert(lap, lap, node, x_periodic, y_periodic, z_periodic);

            break;
        }
        case 3: {
            const int igxn = node->igx;
            const int icxn = node->icx;
            const int igyn = node->igy;
            const int icyn = node->icy;
            const int igzn = node->igz;
            const int iczn = node->icz;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            const int igze = elem->igz;
            const int icze = elem->icz;
            
            const double dx = node->dx;
            const double dy = node->dy;
            const double dz = node->dz;
            
            const double* hcenter  = wcenter;
            
            const double oodx2[3] = {1.0/(dx*dx), 1.0/(dy*dy), 1.0/(dz*dz)};
            const int dis[3][3]   = {{1, icxn, icxn*icyn}, {icxn, icxn*icyn, 1}, {icxn*icyn, 1, icxn}};
            
            const double a00 = 9.0/64.0;
            const double a10 = 3.0/64.0;
            const double a01 = 3.0/64.0;
            const double a11 = 1.0/64.0;
                        
            int i, j, k;
            
            memset(lap, 0.0, node->nc*sizeof(double));
            
            for(k = igze; k < icze - igze; k++) {
                int le   = k * icxe*icye;
                int ln   = k * icxn*icyn;
                
                for(j = igye; j < icye - igye; j++) {
                    int me   = le + j * icxe;
                    int mn   = ln + j * icxn;
                    
                    for(i = igxe; i < icxe - igxe; i++) {
                        int ne = me + i;
                        int nn = mn + i;
                        
                        int nmsw, nmse, nmne, nmnw, npsw, npse, npne, npnw;
                        
                        for (int idim = 0; idim<elem->ndim; idim++) {
                            nmsw = nn;
                            nmse = nn                + dis[idim][1];
                            nmne = nn                + dis[idim][1] + dis[idim][2];
                            nmnw = nn                               + dis[idim][2];
                            npsw = nn + dis[idim][0];
                            npse = nn + dis[idim][0] + dis[idim][1];
                            npne = nn + dis[idim][0] + dis[idim][1] + dis[idim][2];
                            npnw = nn + dis[idim][0]                + dis[idim][2];
                            
                            double dp_sw = (p[npsw] - p[nmsw]);
                            double dp_se = (p[npse] - p[nmse]);
                            double dp_ne = (p[npne] - p[nmne]);
                            double dp_nw = (p[npnw] - p[nmnw]);
                            
                            double flux_sw = hplus[idim][ne] * oodx2[idim] * ( a00 * dp_sw + a10 * dp_se + a11 * dp_ne + a01 * dp_nw);
                            double flux_se = hplus[idim][ne] * oodx2[idim] * ( a01 * dp_sw + a00 * dp_se + a10 * dp_ne + a11 * dp_nw);
                            double flux_ne = hplus[idim][ne] * oodx2[idim] * ( a11 * dp_sw + a01 * dp_se + a00 * dp_ne + a10 * dp_nw);
                            double flux_nw = hplus[idim][ne] * oodx2[idim] * ( a10 * dp_sw + a11 * dp_se + a01 * dp_ne + a00 * dp_nw);
                            
                            lap[nmsw] += flux_sw;
                            lap[nmse] += flux_se;
                            lap[nmne] += flux_ne;
                            lap[nmnw] += flux_nw;
                            lap[npsw] -= flux_sw;
                            lap[npse] -= flux_se;
                            lap[npne] -= flux_ne;
                            lap[npnw] -= flux_nw;
                        }
                        
                        double hc = 0.125 *hcenter[ne];
                        
                        nmsw = nn;
                        nmse = nn             + dis[0][1];
                        nmne = nn             + dis[0][1] + dis[0][2];
                        nmnw = nn                         + dis[0][2];
                        npsw = nn + dis[0][0];
                        npse = nn + dis[0][0] + dis[0][1];
                        npne = nn + dis[0][0] + dis[0][1] + dis[0][2];
                        npnw = nn + dis[0][0]             + dis[0][2];
                        
                        lap[nmsw] += hc*p[nmsw];
                        lap[nmse] += hc*p[nmse];
                        lap[nmne] += hc*p[nmne];
                        lap[nmnw] += hc*p[nmnw];
                        lap[npsw] += hc*p[npsw];
                        lap[npse] += hc*p[npse];
                        lap[npne] += hc*p[npne];
                        lap[npnw] += hc*p[npnw];
                    }
                }
            }
            
#ifdef PRECON
            precon_invert(lap, lap, node, x_periodic, y_periodic, z_periodic);
#endif
            
            if (x_periodic) {
                for(k=igzn; k<iczn-igzn; k++) {
                    int ln = k*icxn*icyn;
                    for(j=igyn; j<icyn-igyn; j++) {
                        int nm = ln + j * icxn + igxn;
                        int np = ln + j * icxn + icxn - igxn - 1;
                        
                        lap[nm]  += lap[np];
                        lap[np]  = 0.0;
                    }
                }
            }
            
            if (y_periodic) {
                for(i=igxn; i<icxn-igxn; i++) {
                    for(k=igzn; k<iczn-igzn; k++) {
                        int ln   = k*icxn*icyn;
                        int nm = ln + i + igyn * icxn;
                        int np = ln + i + (icyn - igyn - 1) * icxn;
                        
                        lap[nm]  += lap[np];
                        lap[np]  = 0.0;
                    }
                }
            }
            
            if (z_periodic) {
                for(j=igyn; j<icyn-igyn; j++) {
                    int mn = j * icxn;
                    for(i=igxn; i<icxn-igxn; i++) {
                        int nm = mn + i + igzn * icxn*icyn;
                        int np = mn + i + (iczn - igzn - 1) * icxn*icyn;
                        
                        lap[nm]  += lap[np];
                        lap[np]  = 0.0;
                    }
                }
            }
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
}

