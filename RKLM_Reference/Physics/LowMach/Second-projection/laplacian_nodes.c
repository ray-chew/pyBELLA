/*******************************************************************************
 File:   laplacian_nodes.c
 Author: Nicola                        Rupert
 Date:   Fri Mar 13 07:56:56 CET 1998  Feb. 2004
 *******************************************************************************/
#include <float.h>
#include <string.h>
#include <math.h>
#include "set_ghostnodes_p.h"
#include "laplacian_nodes.h"
#include "userdata.h"
#include "error.h"
#include "Common.h"
#include "enumerator.h"
#include "ThomasAlgorithmus.h"
#include "mpv.h"
#include "math_own.h"


#ifdef PRECON /* ============================================================================= */

#ifdef PRECON_DIAGONAL_2ND_PROJ 

static enum Boolean diaginv_is_allocated = WRONG;
static double *diaginv;
static double *diag;

#ifdef PRECON_LEGACY
double precon_prepare(
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
    
    const double dx = node->dx;
    const double dy = node->dy;
    const int ndim = node->ndim;
    
    if (ud.precondition == WRONG) {
        for(int nn=0; nn<node->nc; nn++) diag[nn] = diaginv[nn] = 1.0;
        precon_inv_scale = 1.0;
    }
    else {
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
                
                const int icxe = elem->icx;
                const int icye = elem->icy;
                                
                const double* hplusx   = hplus[0];
                const double* hplusy   = hplus[1];
                                                                                 
                const double oodx2 = 1.0 / (dx * dx);
                const double oody2 = 1.0 / (dy * dy);

                int i, j, me, mn, ne, nn, nn1, nnicxn, nn1icxn;
                
                for(nn=0; nn<node->nc; nn++) diag[nn] = diaginv[nn] = 0.0;
                
                for(j = 0; j < icye; j++) {
                    me   = j * icxe;
                    mn   = j * icxn;
                    
                    for(i = 0; i < icxe; i++) {
                        ne       = me + i;
                        
                        nn       = mn + i;
                        nn1      = nn + 1;
                        nnicxn   = nn + icxn;
                        nn1icxn  = nn + 1 + icxn;
                        
                        /* location nn */
                        diag[nn]      += 0.125 * (hplusx[ne]*oodx2+hplusy[ne]*oody2);
                        diag[nn1]     += 0.125 * (hplusx[ne]*oodx2+hplusy[ne]*oody2);
                        diag[nnicxn]  += 0.125 * (hplusx[ne]*oodx2+hplusy[ne]*oody2);
                        diag[nn1icxn] += 0.125 * (hplusx[ne]*oodx2+hplusy[ne]*oody2);
                    }
                }
                
                for (j = igyn; j < icyn-igyn; j++) {mn = j*icxn;
                    for (i = igxn; i < icxn-igxn; i++) {nn = mn + i;
                        diaginv[nn] = 1.0/diag[nn];
                    }
                }
                                
                /* 
                 here we loop over inner nodes only because bdry nodes only pick up to terms
                 in the above scattering loop
                 */
                precon_inv_scale = 0.0;
                for (j = igyn+1; j < icyn-igyn-1; j++) {mn = j*icxn;
                    for (i = igxn+1; i < icxn-igxn-1; i++) {nn = mn + i;
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
    }
    return precon_inv_scale;
}
#else /* PRECON_LEGACY */
double precon_prepare(
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
    
    if (ud.precondition == WRONG) {
        for(int nn=0; nn<node->nc; nn++) diag[nn] = diaginv[nn] = 1.0;
        precon_inv_scale = 1.0;
    }
    else {
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
                
                const double oody  = 0.5 / dy;
                
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
    }
}
#endif /* PRECON_LEGACY */


void precon_apply(
                  double* vec_out,
                  const double* vec_in,
                  const NodeSpaceDiscr *node) {
    for (int nn=0; nn<node->nc; nn++) {
        vec_out[nn] = vec_in[nn]*diag[nn];
    }
}

void precon_invert(
                   double* vec_out,
                   const double* vec_in,
                   const NodeSpaceDiscr *node) {
    for (int nn=0; nn<node->nc; nn++) {
        vec_out[nn] = vec_in[nn]*diaginv[nn];
    }
}
#endif /* PRECON_DIAGONAL_2ND_PROJ =========================================================== */

#ifdef PRECON_VERTICAL_COLUMN_2ND_PROJ /* ==================================================== */
static enum Boolean tridiago_is_allocated = WRONG;
static double *tridiago[3];

double precon_prepare(
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
    
    double precon_inv_scale;
    
    if (tridiago_is_allocated == WRONG) {
        for (int k=0; k<3; k++) {
            tridiago[k] = (double*)malloc(node->nc*sizeof(double));
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
            const int icxn = node->icx;
            
            const int igxe = elem->igx;
            const int icxe = elem->icx;
            const int igye = elem->igy;
            const int icye = elem->icy;
            
            const double dy = node->dy;
            
            const double* hplusy   = hplus[1];
            const double* hcenter  = wcenter;
            
            const double oody  = 0.5 / dy;
            
            const double oody2 = 0.5 / (dy * dy);
            
            double flux_y_left, flux_y_right, hc;
            
            int i, j, me, mn, ne, nn, nn1, nnicxn, nn1icxn;
            
            double pnn      = 1.0;
            double pnn1     = 1.0;
            double pnnicxn  = 1.0;
            double pnn1icxn = 1.0;
            
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
                    
                    /* summand 1.0 takes care of singularity of the tridiagoonal matrix */
                    hc            = 0.25 * (-1.0 + hcenter[ne]);
                                        
                    
                    /* eventually I should transpose the tridiago[][] field for better memory efficiency */
                    /* eventually I should let tridiago[] run from  -1 to +1 */
                    
                    tridiago[1][nn]      += (+ flux_y_left  * ( -pnn      )) + hc*pnn;
                    tridiago[2][nn]      += (+ flux_y_left  * ( +pnnicxn  ));
                    
                    tridiago[1][nn1]     += (+ flux_y_right * ( -pnn1     )) + hc*pnn1;
                    tridiago[2][nn1]     += (+ flux_y_right * ( +pnn1icxn ));
                    
                    tridiago[0][nnicxn]  += (- flux_y_left  * ( -pnn      ));
                    tridiago[1][nnicxn]  += (- flux_y_left  * ( +pnnicxn  )) + hc*pnnicxn;
                    
                    tridiago[0][nn1icxn] += (- flux_y_right * ( -pnn1     ));
                    tridiago[1][nn1icxn] += (- flux_y_right * ( +pnn1icxn )) + hc*pnn1icxn;
                }
            }
            
            This has not been tested !
            precon_inv_scale = 0.0;
            for(nn=0; nn<node->nc; nn++) precon_inv_scale = MAX_own(precon_inv_scale, fabs(1.0/tridiago[1][nn]));
        }
            break;
        case 3: {
            ERROR("function not available");
            break;
        }
            
            
        default: ERROR("ndim not in {1, 2, 3}");
    }

    return precon_inv_scale;
}

void precon_apply(
                  double* vec_out,
                  const double* vec_in,
                  const NodeSpaceDiscr *node) {
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    
    const int igxn = node->igx;
    const int igyn = node->igy;
    
    for (int i = igxn; i<icxn-igxn; i++) {
        for (int j = igyn; j<icyn-igyn; j++) {
            int nn0 = (j-1)*icxn+i;
            int nn1 =   j  *icxn+i;
            int nn2 = (j+1)*icxn+i;
            vec_out[nn1] = tridiago[0][nn1]*vec_in[nn0]+tridiago[1][nn1]*vec_in[nn1]+tridiago[2][nn1]*vec_in[nn2];
        }
    }
}

void precon_invert(
                   double* vec_out,
                   const double* vec_in,
                   const NodeSpaceDiscr *node) {
    
    /* It will be more efficient to immediately store tridiago in the y-first ordering
     to avoid the transposition and intermediate storage in upper, diago, lower,
     but that can be done later
     */
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    
    const int igxn = node->igx;
    const int igyn = node->igy;
    
    const int size = icyn-2*igyn;
    
    double* upper  = (double*)malloc(size*sizeof(double));
    double* diago  = (double*)malloc(size*sizeof(double));
    double* lower  = (double*)malloc(size*sizeof(double));
    double* v_in   = (double*)malloc(size*sizeof(double));
    double* v_out  = (double*)malloc(size*sizeof(double));
    
    for (int i=igxn; i<icxn-igxn; i++) {
        for (int j=igyn; j<icyn-igyn; j++) {
            int j_inn = j-igyn;
            int nn    = j*icxn+i;
            lower[j_inn] = tridiago[0][nn];
            diago[j_inn] = tridiago[1][nn];
            upper[j_inn] = tridiago[2][nn];
            v_in[j_inn]  = vec_in[nn];
        }
        Thomas_Algorithm(v_out, v_in, upper, diago, lower, size);
        for (int j=igyn; j<icyn-igyn; j++) {
            int j_inn = j-igyn;
            int nn    = j*icxn+i;
            vec_out[nn] = v_out[j_inn];
        }
    }

    free(upper);
    free(diago);
    free(lower);
    free(v_in );
    free(v_out);
    
}
#endif /* PRECON_VERTICAL_COLUMN_2ND_PROJ ==================================================== */

#else  /* PRECON ============================================================================= */
double precon_prepare(
                           const NodeSpaceDiscr* node,
                           const ElemSpaceDiscr* elem,
                           const double* hplus[3],
                           const double* wcenter,
                           const int x_periodic,
                           const int y_periodic,
                           const int z_periodic) {
}

void precon_apply(
                  double* vec_out,
                  const double* vec_in,
                  const NodeSpaceDiscr *node) {
}

void precon_invert(
                   double* vec_out,
                   const double* vec_in,
                   const NodeSpaceDiscr *node) {
}
#endif /* PRECON ============================================================================= */

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
            
            precon_invert(lap, lap, node);
            
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
            
            assert(0); /* implicit gravity not implemented for 3D yet */
            
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
            precon_invert(lap, lap, node);
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

