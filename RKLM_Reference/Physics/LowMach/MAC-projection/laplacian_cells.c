#/*******************************************************************************
File:   laplacian_cells.c
Author: Nicola                        Rupert     
Date:   Fri Mar 13 07:56:56 CET 1998  Feb. 2004
*******************************************************************************/

#include <assert.h> 
#include <math.h>
#include <float.h>
#include <string.h>

#include "Common.h"
#include "set_ghostcells_p.h"
#include "laplacian_cells.h"
#include "error.h"
#include "userdata.h"
#include "thermodynamic.h"
#include "EOS.h"
#include "math_own.h"

/* ========================================================================== */
#ifdef PRECON

#ifdef PRECON_DIAGONAL_1ST_PROJ
static enum Boolean diaginv_c_is_allocated = WRONG;
static double *diaginv_c;
static double *diag_c;

void precon_c_prepare(
                      const NodeSpaceDiscr* node,
                      const ElemSpaceDiscr* elem,
                      const double* hplus[3],
                      const double* hcenter,
                      const double* hgrav) {
    
    /*
     the following defines diag/diaginv directly through the calculations in
     EnthalpyWeightedLap_bilinear_p()
     As a consequence, also the boundary elements have the real diagonal values,
     as opposed to what I did previously, when all the elements received the
     diagonal value from a full stencil.
     */
    
    
    if (diaginv_c_is_allocated == WRONG) {
        diag_c    = (double*)malloc(node->nc*sizeof(double));
        diaginv_c = (double*)malloc(node->nc*sizeof(double));
        diaginv_c_is_allocated = CORRECT;
    }
    
    const int ndim = elem->ndim;
    double diagscale;
    
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
            const double oody   = 1.0 / dy;
            const double oodx2  = 1.0 / (dx * dx);
            const double oody2  = 1.0 / (dy * dy);
            
            const double pnc    = 1.0;
            
            const double* hplusx = hplus[0];
            const double* hplusy = hplus[1];
            const double* hc     = hcenter;
            const double* hg     = hgrav;
            
            for(int nc=0; nc<elem->nc; nc++) diag_c[nc] = diaginv_c[nc] = 0.0;
            
            for(int j = igy; j < icy - igy; j++) {int mc = j * icx;
                for(int i = igx; i < icx - igx; i++) {int nc  = mc + i;
                    
                    int ox  = j * ifx + i;
                    int oy  = i * ify + j;
                    int o_e = ox + 1;
                    int o_w = ox;
                    int o_n = oy + 1;
                    int o_s = oy;
                    
                    diag_c[nc]  = oodx2 * ( (  hplusx[o_e] * (-pnc) - hplusx[o_w] * pnc
                                         )
                                       + 0.125 * (  hplusx[o_e] * (   pnc )
                                                  - hplusx[o_w] * ( - pnc )
                                                  - hplusx[o_e] * ( - pnc )
                                                  + hplusx[o_w] * (   pnc )
                                                  )
                                       );
                    
                    diag_c[nc] += oody2 * ( (  hplusy[o_n] * ( - pnc ) - hplusy[o_s] * ( pnc )
                                         )
                                       + 0.125 * (  hplusy[o_n] * (  pnc )
                                                  - hplusy[o_s] * (- pnc )
                                                  - hplusy[o_n] * (- pnc )
                                                  + hplusy[o_s] * (  pnc )
                                                  )
                                       );
                    
                    diag_c[nc] += hc[nc] * pnc;
                    
                    /* implicit gravity works only for gravity in y-direction */
                    diag_c[nc] += oody * (hg[o_n] * 0.5*( pnc ) - hg[o_s] * 0.5*( pnc ));
                }
            }

            diagscale = 3.0 /(elem->dx*elem->dx) + 3.0 /(elem->dy*elem->dy);

            for (int j = igy; j < icy-igy; j++) {int mc = j*icx;
                for (int i = igx; i < icx-igx; i++) {int nc = mc + i;
                    diag_c[nc]   *= diagscale;
                    diaginv_c[nc] = 1.0/diag_c[nc];
                }
            }

            break;
        }
        case 3: {
            ERROR("function not available for 3D");
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }

}

void precon_c_apply(
                    double* vec_out,
                    const double* vec_in,
                    const ElemSpaceDiscr *elem) {
    for (int nc=0; nc<elem->nc; nc++) {
        vec_out[nc] = vec_in[nc]*diag_c[nc];
    }
}

void precon_c_invert(
                     double* vec_out,
                     const double* vec_in,
                     const ElemSpaceDiscr *elem) {
    for (int nc=0; nc<elem->nc; nc++) {
        vec_out[nc] = vec_in[nc]*diaginv_c[nc];
    }
}
#endif /* PRECON_DIAGONAL */

#ifdef PRECON_VERTICAL_COLUMN_1ST_PROJ
Nothing defined yet
#endif /* PRECON_VERTICAL_COLUMN */

#else  /* PRECON */
void precon_c_prepare(
                      const NodeSpaceDiscr* node,
                      const ElemSpaceDiscr* elem,
                      const double* hplus[3],
                      const double* hcenter,
                      const double* hgrav) {
}

void precon_c_apply(
                    double* vec_out,
                    const double* vec_in,
                    const ElemSpaceDiscr *elem) {
}

void precon_c_invert(
                     double* vec_out,
                     const double* vec_in,
                     const ElemSpaceDiscr *elem) {
}
#endif /* PRECON */

 
#ifdef SHIFTED_COEFFICIENTS_PROJ1
/* ========================================================================== */

void EnthalpyWeightedLap_bilinear_p(
									const ElemSpaceDiscr* elem, 
                                    const NodeSpaceDiscr* node,
									const double* p,
									const double* hplus[3],
									const double* hcenter,
									const double* hgrav,
									const ConsVars* Sol,
									const MPV* mpv, 
									const double dt,
									const double theta,
									double* lap) {
	
	const int ndim = elem->ndim;
	
	switch(ndim) {
		case 1: {    
			ERROR("function not available");
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

			/* for MG-scaling with elem->scale_factor; see old version of routine */

            const double oody      = 1.0 / dy;

            const double oodx2     = 1.0 / (dx * dx);
			const double oody2     = 1.0 / (dy * dy);
            
			const double* hplusx   = hplus[0];
			const double* hplusy   = hplus[1];
			const double* hc       = hcenter;
			const double* hg       = hgrav;
            
            double dlapsum;
            
            dlapsum = 0.0;
			
			int n_c, n_e, n_ne, n_n, n_nw, n_w, n_sw, n_s, n_se;
			int o_en, o_ec, o_es, o_nw, o_nc, o_ne, o_wn, o_wc, o_ws, o_se, o_sc, o_sw, ox, oy;
			
			int i, j, m, n;
            			
			for(j = igy; j < icy - igy; j++) {
                m = j * icx;				
				
                for(i = igx; i < icx - igx; i++) {
                    n  = m + i;
                    
					ox = j * ifx + i;
					oy = i * ify + j;
					
					n_c    = n;
					n_e    = n + 1      ;
					n_ne   = n + 1 + icx;
					n_n    = n     + icx;
					n_nw   = n - 1 + icx;
					n_w    = n - 1      ;
					n_sw   = n - 1 - icx;
					n_s    = n     - icx;
					n_se   = n + 1 - icx;
					
                    o_es   = ox + 1 - ifx;
                    o_ec   = ox + 1;
					o_en   = ox + 1 + ifx;

                    o_ws   = ox - ifx;
                    o_wc   = ox;
                    o_wn   = ox + ifx;
					
                    o_ne   = oy + 1 + ify;
                    o_nc   = oy + 1;
					o_nw   = oy + 1 - ify;
                    
                    o_se   = oy + ify;
                    o_sc   = oy;
					o_sw   = oy - ify;
					
					lap[n]  = oodx2 * (  0.75  * (  hplusx[o_ec] * (p[n_e ] - p[n_c ]) - hplusx[o_wc] * (p[n_c ] - p[n_w ])) 
                                       + 0.125 * (  hplusx[o_es] * (p[n_se] - p[n_s ]) - hplusx[o_ws] * (p[n_s ] - p[n_sw])  
                                                  + hplusx[o_en] * (p[n_ne] - p[n_n ]) - hplusx[o_wn] * (p[n_n ] - p[n_nw])   
                                                 )
                                       );
                    
					lap[n] += oody2 * (  0.75  * (  hplusy[o_nc] * (p[n_n ] - p[n_c ]) - hplusy[o_sc] * (p[n_c ] - p[n_s ]))
                                       + 0.125 * (  hplusy[o_ne] * (p[n_ne] - p[n_e ]) - hplusy[o_se] * (p[n_e ] - p[n_se]) 
                                                  + hplusy[o_nw] * (p[n_nw] - p[n_w ]) - hplusy[o_sw] * (p[n_w ] - p[n_sw]) 
                                                  )
                                       );
                    
                    lap[n] += oody * 0.5 * (hg[o_nc] * (p[n_n] + p[n_c]) - hg[o_sc] * (p[n_s] + p[n_c]));
                    
                    /*
                     {
                     double dlap1, dlap2;
                     
                     dlap1 = oody * 0.5 * (hg[o_nc] * (p[n_n] + p[n_c]) - hg[o_sc] * (p[n_s] + p[n_c]));
                    
                     dlap2 = oody * 0.5 * ( 0.75  * (hg[o_nc] * (p[n_n ] + p[n_c ]) - hg[o_sc] * (p[n_s ] + p[n_c ]))
                                          +0.125 * (hg[o_ne] * (p[n_ne] + p[n_e ]) - hg[o_se] * (p[n_se] + p[n_e ]))
                                          +0.125 * (hg[o_nw] * (p[n_nw] + p[n_w ]) - hg[o_sw] * (p[n_sw] + p[n_w ]))
                                                    );
                     dlapsum += fabs(dlap1 - dlap2);

                     lap[n] += dlap1;
                     }
                    */
 					
					lap[n] += hc[n] * p[n_c];
                    
				}
			}    
            
            /* printf("dlap1 = %e\t dlap2 = %e\t dlapsum = %e\n", dlap1, dlap2, dlapsum); */
			
#ifdef PRECON
            precon_c_invert(lap, lap, elem);
#endif
			break;
		}
        case 3:  ERROR("SHIFTED_COEFFICIENTS Version not implemented in 3D yet");
		default: ERROR("ndim not in {1, 2, 3}");
	}
}

#else /* SHIFTED_COEFFICIENTS */

/* ========================================================================== */

void EnthalpyWeightedLap_bilinear_p(
                                    const ElemSpaceDiscr* elem, 
                                    const NodeSpaceDiscr* node,
                                    const double* p,
                                    const double* hplus[3],
                                    const double* hcenter,
                                    const double* hgrav,
                                    const ConsVars* Sol,
                                    const MPV* mpv, 
                                    const double dt,
                                    const double theta,
                                    double* lap) {
    
    const int ndim = elem->ndim;
    
    switch(ndim) {
        case 1: {    
            ERROR("function not available");
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
            
            /* for MG-scaling with elem->scale_factor; see old version of routine */
            
            const double oody      = 1.0 / dy;
            
            const double oodx2     = 1.0 / (dx * dx);
            const double oody2     = 1.0 / (dy * dy);
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            const double* hc       = hcenter;
            const double* hg       = hgrav;
            
            int n_c, n_e, n_ne, n_n, n_nw, n_w, n_sw, n_s, n_se;
            int o_e, o_n, o_w, o_s, ox, oy;
            
            int i, j, m, n;
            
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
                    
                    lap[n]  = oodx2 * ( (  hplusx[o_e] * (p[n_e ] - p[n_c ]) - hplusx[o_w] * (p[n_c ] - p[n_w ]) 
                                         ) 
                                       + 0.125 * (  hplusx[o_e] * (  (p[n_ne] - p[n_n ] ) - (p[n_e ] - p[n_c ] ) ) 
                                                  - hplusx[o_w] * (  (p[n_n ] - p[n_nw] ) - (p[n_c ] - p[n_w ] ) )  
                                                  - hplusx[o_e] * (  (p[n_e ] - p[n_c ] ) - (p[n_se] - p[n_s ] ) ) 
                                                  + hplusx[o_w] * (  (p[n_c ] - p[n_w ] ) - (p[n_s ] - p[n_sw] ) ) 
                                                  )
                                       );
                    
                    lap[n] += oody2 * ( (  hplusy[o_n] * ( p[n_n ] - p[n_c ] ) - hplusy[o_s] * ( p[n_c ] - p[n_s ] ) 
                                         )
                                       + 0.125 * (  hplusy[o_n] * (  (p[n_ne] - p[n_e ] ) - (p[n_n ] - p[n_c ] ) ) 
                                                  - hplusy[o_s] * (  (p[n_e ] - p[n_se] ) - (p[n_c ] - p[n_s ] ) ) 
                                                  - hplusy[o_n] * (  (p[n_n ] - p[n_c ] ) - (p[n_nw] - p[n_w ] ) ) 
                                                  + hplusy[o_s] * (  (p[n_c ] - p[n_s ] ) - (p[n_w ] - p[n_sw] ) ) 
                                                  )
                                       );
                    lap[n] += oody * (hg[o_n] * 0.5*(p[n_n] + p[n_c]) - hg[o_s] * 0.5*(p[n_s] + p[n_c]));
                    
                    lap[n] += hc[n] * p[n_c];
                    
                }
            }    
            
            
#ifdef PRECON
            precon_c_invert(lap, lap, elem);
#endif
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
            
#ifdef PRECON
            precon_c_invert(lap, lap, elem);
#endif
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
}

#endif /* SHIFTED_COEFFICIENTS */
