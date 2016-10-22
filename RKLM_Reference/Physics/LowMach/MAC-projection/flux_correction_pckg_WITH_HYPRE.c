/*******************************************************************************
 File:   flux_correction.c
 Author: Rupert
 Date:   Feb  14  2004
 *******************************************************************************/
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "Common.h"
#include "set_ghostcells_p.h"
#include "mpv.h"
#include "memory.h"
#include "flux_correction.h"
#include "variable_coefficient_poisson_cells.h"
#include "BICGSTAB.h"
#include "variable.h"
#include "warning.h"
#include "enumerator.h"
#include "error.h"
#include "io.h"
#include "userdata.h"
#include "Eos.h"
#include "thermodynamic.h"
#include "math_own.h"
#include "recovery.h"
#include "laplacian_cells.h"

static enum Constraint integral_condition(
										  ConsVars* flux[3],
										  double* rhs, 
										  ConsVars* Sol,
										  const double dt,
										  const ElemSpaceDiscr* elem,
										  MPV* mpv);

static void controlled_variable_change_explicit(
												double* rhs, 
												const ElemSpaceDiscr* elem, 
												ConsVars* Sol_new, 
												ConsVars* Sol_old, 
												double dt,
												const MPV* mpv);

static void rhs_fix_for_open_boundaries(
										double* rhs, 
										const ElemSpaceDiscr* elem, 
										ConsVars* Sol_new, 
										ConsVars* Sol_old, 
										ConsVars* flux[3],
										double dt,
										MPV* mpv);

static void flux_fix_for_open_boundaries(
										 ConsVars* flux[3],
										 const ElemSpaceDiscr* elem, 
										 MPV* mpv);

static void flux_correction_due_to_pressure_gradients(
													  ConsVars* flux[3],
													  VectorField* buoy,
													  const ElemSpaceDiscr* elem,
													  ConsVars* Sol,
													  ConsVars* Sol0,
													  const MPV* mpv, 
													  double* hplus[3],
													  double* hgrav,
													  double* dp2,
													  const double t,
													  const double dt,
													  const double implicitness);

static void flux_correction_due_to_pressure_values(
												   ConsVars* flux[3],
												   VectorField* buoy,
												   const ElemSpaceDiscr* elem,
												   ConsVars* Sol,
												   double* dp2, 
												   const double dt);

/* ========================================================================== */

#if 1
static int rhs_output_count = 0;
#endif

void flux_correction(
					 ConsVars* flux[3],
					 VectorField* buoy,
					 const ElemSpaceDiscr* elem_base_grid,
					 ConsVars* Sol, 
					 ConsVars* Sol0, 
					 const double t,
					 const double dt,
					 const double implicitness,
                     const int step) {
	
    extern User_Data ud;
	extern MPV* mpv;
	
	const ElemSpaceDiscr* elem = mpv->Level[0]->elem;
    const NodeSpaceDiscr* node = mpv->Level[0]->node;

	double** hplus       = mpv->Level[0]->wplus;
	double*  hcenter     = mpv->Level[0]->wcenter;
	double*  hgrav       = mpv->Level[0]->wgrav;
	
	double* rhs          = mpv->Level[0]->rhs;
	double* dp2		     = mpv->Level[0]->p;
    
    rhs_output_count     = 0;
    
	int n;
        
    printf("\n\n====================================================");
    printf("\nFirst Projection");
    printf("\n====================================================\n");
	    
	operator_coefficients(hplus, hcenter, hgrav, elem, Sol, Sol0, mpv, dt);
    
	/* obtain finest grid data for the MG hierarchy */
	controlled_variable_change_explicit(rhs, elem, Sol, Sol0, dt, mpv);
    
    assert(integral_condition(flux, rhs, Sol, dt, elem, mpv) != VIOLATED); 
    rhs_fix_for_open_boundaries(rhs, elem, Sol, Sol0, flux, dt, mpv);
    
#if 1
    extern User_Data ud;
    FILE *prhsfile = NULL;
    char fn2[200], fieldname2[90];
    
    sprintf(fn2, "%s/rhs_cells/rhs_cells_00%d.hdf", ud.file_name, rhs_output_count);
    rhs_output_count++;
    sprintf(fieldname2, "rhs_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             rhs,
             fn2,
             fieldname2);
    
    /*
    sprintf(fn2, "%s/Tests/p2_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "p2_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             dp2,
             fn2,
             fieldname2);

    sprintf(fn2, "%s/Tests/hx_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hx_c");
    
    WriteHDF(prhsfile,
             2*(mpv->Level[0]->elem->icy-2*mpv->Level[0]->elem->igy),
             mpv->Level[0]->elem->ifx-2*mpv->Level[0]->elem->igx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             h[0],
             fn2,
             fieldname2);
    
    sprintf(fn2, "%s/Tests/hy_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hy_c");
    
    WriteHDF(prhsfile,
             2*(mpv->Level[0]->elem->ify-2*mpv->Level[0]->elem->igy),
             mpv->Level[0]->elem->icx-2*mpv->Level[0]->elem->igx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             h[1],
             fn2,
             fieldname2);
    
    sprintf(fn2, "%s/Tests/hplusx_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hplusx_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->ifx,
             mpv->Level[0]->elem->icy,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             hplus[0],
             fn2,
             fieldname2);
    
    sprintf(fn2, "%s/Tests/hplusy_c_test.hdf", ud.file_name);
    sprintf(fieldname2, "hplusy_c");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->ify,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             hplus[1],
             fn2,
             fieldname2);
    
 */
#endif

    
    /**/
    variable_coefficient_poisson_cells(dp2, rhs, (const double **)hplus, hcenter, hgrav, Sol, elem, node, implicitness);

    set_ghostcells_p2(dp2, (const double **)hplus, hgrav, elem, elem->igx);

    flux_correction_due_to_pressure_gradients(flux, buoy, elem, Sol, Sol0, mpv, hplus, hgrav, dp2, t, dt, implicitness);
    if (ud.p_flux_correction) {
        flux_correction_due_to_pressure_values(flux, buoy, elem, Sol, dp2, dt); 
    }

    flux_fix_for_open_boundaries(flux, elem, mpv);  

    memcpy(mpv->dp2_cells, dp2, elem->nc*sizeof(double));

    /* store results in mpv-fields */        
    for(n=0; n<elem->nc; n++) mpv->p2_cells[n] = mpv->p2_cells[n] + dp2[n];
    
	set_ghostcells_p2(mpv->p2_cells, (const double **)hplus, hgrav, elem, elem->igx);
}

/* ========================================================================== */
 
static void controlled_variable_change_explicit(
                                                double* rhs, 
                                                const ElemSpaceDiscr* elem, 
                                                ConsVars* Sol_new, 
                                                ConsVars* Sol_old, 
                                                double dt,
                                                const MPV* mpv) {
        
    memset(rhs, 0.0, elem->nc*sizeof(double));
    
    const int igx = elem->igx;
    const int icx = elem->icx;
    const int igy = elem->igy;
    const int icy = elem->icy;
    const int igz = elem->igz;
    const int icz = elem->icz;
    
    const double factor = 2.0 / (dt * dt);
        
    int i, j, k, l, m, n;
    
    for(i=0; i<elem->nc; i++) rhs[i] = 0.0;
        
    for(k = igz; k < icz - igz; k++) {l = k * icx*icy;
        for(j = igy; j < icy - igy; j++) {m = l + j * icx;
            for(i = igx; i < icx - igx; i++) {n = m + i;
                double rhs_value = - factor * (Sol_new->rhoY[n] - Sol_old->rhoY[n]);
                rhs[n]   = rhs_value;
            }
        }
    }
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_cells.hdf");
    sprintf(fieldname, "rhs-cells");    
    WriteHDF(prhsfile, elem->icx, elem->icy, elem->icz, elem->ndim, rhs, fn, fieldname);
#endif
}

/* ========================================================================== */

static enum Constraint integral_condition(
										  ConsVars* flux[3],
										  double* rhs,
										  ConsVars* Sol,
										  const double dt,
										  const ElemSpaceDiscr* elem,
										  MPV* mpv) {
	
	const int ndim = elem->ndim;
	
	switch(ndim) {
			
		case 1: {
			ERROR("function not available");
			break;
		}
			
		case 2: {
			int i, j, jifx, iify;
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int ifx = elem->ifx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int ify = elem->ify;
			const ConsVars* f = flux[0];
			const ConsVars* g = flux[1];
			double F, G;
			double F_rho, G_rho;
			double F_rhoY, G_rhoY;
						
			/* this implementation does not allow for non-zero global flux divergence */
			for(j = igy, F = 0.0, F_rho = 0.0, F_rhoY = 0.0; j < icy - igy; j++) {
				jifx = j * ifx;
				F      += f->rhoe[jifx + icx - igx] - f->rhoe[jifx +  igx];
				F_rho  += f->rho[jifx + icx - igx]  - f->rho[jifx +  igx];
				F_rhoY += f->rhoY[jifx + icx - igx] - f->rhoY[jifx +  igx];
			} 
			for(i = igx, G = 0.0, G_rho = 0.0, G_rhoY = 0.0; i < icx - igx; i++) {
				iify =  i * ify;
				G      += g->rhoe[iify + icy - igy] - g->rhoe[iify + igy];
				G_rho  += g->rho[iify + icy - igy]  - g->rho[iify + igy];
				G_rhoY += g->rhoY[iify + icy - igy] - g->rhoY[iify + igy];
			}  
						
            /*
            if(0){
				printf("F_rhoY = %e, G_rhoY = %e, DBL_EPSILON = %e\n", F_rhoY, G_rhoY, DBL_EPSILON);
				return VIOLATED;
			}
			else
             */
            {
				double tmp = 0.0;
				double rhs_max = 0.0;
				double threshold = 1000*DBL_EPSILON;
				int cnt = 0;
				int i,j,m,n;
				for(j = igy; j < icy - igy; j++) {m = j*icx;
					for(i = igx; i < icx - igx; i++) {n = m + i;
						tmp += rhs[n];
						rhs_max = MAX_own(fabs(rhs[n]),rhs_max);
						cnt++;
					}
				}
				tmp /= cnt;
				if(fabs(tmp) < threshold) {
					printf("integral_condition_cells = OK; tmp = %e, \n F = %e, F_rho = %e, F_rhoY = %e, \n G = %e, G_rho = %e, G_rhoY = %e\n", fabs(tmp), F, F_rho, F_rhoY, G, G_rho, G_rhoY);
					return(SATISFIED);
				}
				else {
					printf("integral_condition_cells = VIOLATED; tmp = %e, rhs_max = %e\n F = %e, F_rho = %e, F_rhoY = %e \n G = %e, G_rho = %e, G_rhoY = %e\n", fabs(tmp), rhs_max, F, F_rho, F_rhoY, G, G_rho, G_rhoY);
					return(SATISFIED);
			
                }
			}	
			break;
		}
			
		case 3: {
			int i, j, k, l, m, n;
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int ifx = elem->ifx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int ify = elem->ify;
			const int igz = elem->igz;
			const int icz = elem->icz;
			const int ifz = elem->ifz;
			const int ifxicy = ifx * icy;
			const int ifyicz = ify * icz;
			const int ifzicx = ifz * icx;
			const ConsVars* f = flux[0];
			const ConsVars* g = flux[1];
			const ConsVars* h = flux[2];
			double F, G, H;
			
			printf("\nSanity check for rho-flux and rhoe-flux not implemented in 3D \n");
			
			for(k = igz, F = 0.0; k < icz - igz; k++) {l = k * ifxicy;
				for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
					F += f->rhoY[m + icx - igx] - f->rhoY[m +  igx];
				} 
			}
			for(i = igx, G = 0.0; i < icx - igx; i++) {m = i * ifyicz;
				for(k = igz; k < icz - igz; k++) {n = m + k * ify;
					G += g->rhoY[n + icy - igy] - g->rhoY[n + igy];
				} 
			}
			for(j = igy, H = 0.0; j < icy - igy; j++) {n = j * ifzicx;
				for(i = igx; i < icx - igx; i++) {l = n + i * ifz;
					H += h->rhoY[l + icz - igz] - h->rhoY[l + igz];
				}
			}  
			if(fabs(F * elem->dy + G * elem->dx * H * elem->dz) > 100*DBL_EPSILON) 
				return VIOLATED;
			else
				return SATISFIED;
			break;
		}
			
		default: ERROR("ndim not in {1, 2, 3}");
			return(VIOLATED);
			
    }
	
	return(VIOLATED);
}

/* ========================================================================== */

void operator_coefficients(
                           double* hplus[3], 
                           double* wcenter,
                           double* wgrav,
                           const ElemSpaceDiscr* elem,
                           const ConsVars* Sol,
                           const ConsVars* Sol0,
                           const MPV* mpv, 
                           const double dt) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	
    const double Gammainv = th.Gammainv;
    
	const int ndim = elem->ndim;
	
    const int impl_grav_th = ud.implicit_gravity_theta;
    const int impl_grav_pr = ud.implicit_gravity_press;
    
	const double implicitness = ud.implicitness;
    const double ccw = 2.0; /* ccenterweight   4.0 */
	const double ccenter = - ccw * (ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt); 
    /* for ccw = 4.0: 
       first factor 2.0 from reinterpretation of Helmholtz-solution as 
       (p^{n+1}-p^{n}) and SECOND-ORDER update to t^{n+1/2} in the first projection 
     */
	const double cexp    = 1.0-th.gamm;
	
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
			double* hx = hplus[0];
			double* hy = hplus[1];
			double* hc = wcenter;
			double* hg = wgrav;
			
			double hi, him, hj, hjm, g, gimp, thet, thetm, Msq;
			
			int i, j, m, n, ic, icm, jc, jcm;
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[0];
            
			for(j = igy-1; j < icy - igy+1; j++) {
                m = j * ifx;
				
                for(i = igx-1; i < ifx - igx+1; i++) {
                    n     = m + i;
                    ic    = n - j;
					icm   = ic - 1; 
					
#ifdef THERMCON
                    hi    = 0.5 * ( Sol->rhoY[ic] * Sol->rhoY[ic] / Sol->rho[ic] + Sol0->rhoY[ic] * Sol0->rhoY[ic] / Sol0->rho[ic] ) * Gammainv;   
                    him   = 0.5 * ( Sol->rhoY[icm] * Sol->rhoY[icm] / Sol->rho[icm] + Sol0->rhoY[icm] * Sol0->rhoY[icm] / Sol0->rho[icm] ) * Gammainv;
#else
                    hi    = 0.5 * ( Sol->rhoY[ic] / Sol->rho[ic] + Sol0->rhoY[ic] / Sol0->rho[ic] );   
					him   = 0.5 * ( Sol->rhoY[icm] / Sol->rho[icm] + Sol0->rhoY[icm] / Sol0->rho[icm] );
#endif
                    thet  = Sol->rhoY[ic]  / Sol->rho[ic] ;
                    thetm = Sol->rhoY[icm] / Sol->rho[icm];
                    
                    gimp  = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dx*0.5*(thet+thetm)));
                    
					hx[n] = 0.5 * (hi + him) * implicitness * gimp;
                    
					assert(hx[n] > 0.0);
				}
			}
			
            
            g = ud.gravity_strength[1];
            
			for(i = 0; i < icx; i++) {
                n = i * ify;
				
                for(j = igy-1; j < ify - igy+1; j++) {
                    m     = n + j;
                    jc    = j * icx + i;
					jcm   = jc - icx;          
					
#ifdef THERMCON
                    hj    = 0.5 * ( Sol->rhoY[jc] * Sol->rhoY[jc] / Sol->rho[jc] + Sol0->rhoY[jc] * Sol0->rhoY[jc] / Sol0->rho[jc]) * Gammainv;
                    hjm   = 0.5 * ( Sol->rhoY[jcm] * Sol->rhoY[jcm] / Sol->rho[jcm] + Sol0->rhoY[jcm] * Sol0->rhoY[jcm] / Sol0->rho[jcm]) * Gammainv;
#else
					hj    = 0.5 * ( Sol->rhoY[jc] / Sol->rho[jc] + Sol0->rhoY[jc] / Sol0->rho[jc]);
					hjm   = 0.5 * ( Sol->rhoY[jcm] / Sol->rho[jcm] + Sol0->rhoY[jcm] / Sol0->rho[jcm]);
#endif
                    
                    thet  = Sol->rhoY[jc]  / Sol->rho[jc] ;
                    thetm = Sol->rhoY[jcm] / Sol->rho[jcm];
                     
                    gimp  = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dy*0.5*(thet+thetm)));
                    
					hy[m] = 0.5 * (hj + hjm) * implicitness * gimp;
                    
                    hg[m] = th.gamminv * pow(0.5*(Sol->rhoY[jc]+Sol->rhoY[jcm]),cexp) * g * impl_grav_pr * gimp; 
                    /* hg[m] = th.gamminv * pow(0.5*(Sol0->rhoY[jc]+Sol0->rhoY[jcm]),cexp) * (g/Msq) * impl_grav_pr * gimp; */
                    
					assert(hy[m] > 0.0);
				}
			}
            
			for(j = igy; j < icy - igy; j++) {m = j * icx;
				for(i = igx; i < icx - igx; i++) {n = m + i;
					hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
				}
			}
			
			
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
			double* hx = hplus[0];
			double* hy = hplus[1];
			double* hz = hplus[2];
			double* hc = wcenter;
			double* hg = wgrav;
			
			double hi, him, hj, hjm, hk, hkm, g, gimp, thet, thetm, Msq;
			
			int i, j, k, l, m, n, ic, icm, jc, jcm, kc, kcm;
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[0];
            assert(g==0.0); /* implicit gravity only for y-direction */
            
			for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                    for(i = igx; i < ifx - igx; i++) {n = m + i;
                        ic  = k*icx*icy + j*icx + i;
                        icm = ic - 1; 
#ifdef THERMCON
                        hi    = 0.5 * (Sol->rhoY[ic] *Sol->rhoY[ic] /Sol->rho[ic]  + Sol0->rhoY[ic] *Sol0->rhoY[ic] /Sol0->rho[ic] ) * Gammainv;   
                        him   = 0.5 * (Sol->rhoY[icm]*Sol->rhoY[icm]/Sol->rho[icm] + Sol0->rhoY[icm]*Sol0->rhoY[icm]/Sol0->rho[icm]) * Gammainv;
#else
                        hi    = 0.5 * (Sol->rhoY[ic] /Sol->rho[ic]  + Sol0->rhoY[ic] /Sol0->rho[ic] );   
                        him   = 0.5 * (Sol->rhoY[icm]/Sol->rho[icm] + Sol0->rhoY[icm]/Sol0->rho[icm]);
#endif
					
                        thet  = Sol->rhoY[ic]  / Sol->rho[ic] ;
                        thetm = Sol->rhoY[icm] / Sol->rho[icm];
                        gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dx*0.5*(thet+thetm)));
                    
                        hx[n] = 0.5 * (hi + him) * implicitness * gimp;
                        assert(hx[n] > 0.0);
                    }
                }
            }
			
            
            g = ud.gravity_strength[1];
            
			for(i = igx; i < icx - igx; i++) {l = i * ify*icz;
                for(k = igz; k < icz - igz; k++) {m = l + k * ify;
                    for(j = igy; j < ify - igy; j++) {n = m + j;
                        jc  = k*icx*icy + j*icx + i;
                        jcm = jc - icx;          
#ifdef THERMCON
                        hj       = 0.5 * (Sol->rhoY[jc] *Sol->rhoY[jc] /Sol->rho[jc] + Sol0->rhoY[jc] *Sol0->rhoY[jc] /Sol0->rho[jc] ) * Gammainv;
                        hjm      = 0.5 * (Sol->rhoY[jcm]*Sol->rhoY[jcm]/Sol->rho[jcm]+ Sol0->rhoY[jcm]*Sol0->rhoY[jcm]/Sol0->rho[jcm]) * Gammainv;
#else
                        hj       = 0.5 * (Sol->rhoY[jc] /Sol->rho[jc] + Sol0->rhoY[jc] /Sol0->rho[jc] );
                        hjm      = 0.5 * (Sol->rhoY[jcm]/Sol->rho[jcm]+ Sol0->rhoY[jcm]/Sol0->rho[jcm]);
#endif
					
                        thet  = Sol->rhoY[jc]  / Sol->rho[jc] ;
                        thetm = Sol->rhoY[jcm] / Sol->rho[jcm];
                        gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dy*0.5*(thet+thetm)));
                    
                        hy[n] = 0.5 * (hj + hjm) * implicitness * gimp;
                        
                        hg[m] = th.gamminv * pow(0.5*(Sol->rhoY[jc]+Sol->rhoY[jcm]),cexp) * g * implicitness * gimp;
                        /* hg[m] = th.gamminv * pow(0.5*(Sol->rhoY[jc]+Sol->rhoY[jcm]),cexp) * (g/Msq) * implicitness * gimp; */

                        assert(hy[n] > 0.0);
                    }
                }
            }

            g = ud.gravity_strength[2];
            assert(g==0.0); /* implicit gravity only for y-direction */
            
			for(j = igy; j < icy - igy; j++) {l = j * ifz*icx;
				for(i = igx; i < icx - igx; i++) {m = l + i * ifz;
                    for(k = igz; k < ifz - igz; k++) {n = m + k;
                        kc  = k*icx*icy + j*icx + i;
                        kcm = kc - icx*icy;          
#ifdef THERMCON
                        hk       = 0.5 * (Sol->rhoY[kc] *Sol->rhoY[kc] /Sol->rho[kc] + Sol0->rhoY[kc] *Sol0->rhoY[kc] /Sol0->rho[kc] ) * Gammainv;
                        hkm      = 0.5 * (Sol->rhoY[kcm]*Sol->rhoY[kcm]/Sol->rho[kcm]+ Sol0->rhoY[kcm]*Sol0->rhoY[kcm]/Sol0->rho[kcm]) * Gammainv;
#else
                        hk       = 0.5 * (Sol->rhoY[kc] /Sol->rho[kc] + Sol0->rhoY[kc] /Sol0->rho[kc] );
                        hkm      = 0.5 * (Sol->rhoY[kcm]/Sol->rho[kcm]+ Sol0->rhoY[kcm]/Sol0->rho[kcm]);
#endif
					
                        thet  = Sol->rhoY[kc]  / Sol->rho[kc] ;
                        thetm = Sol->rhoY[kcm] / Sol->rho[kcm];
                        gimp   = 1.0 / (1.0 + impl_grav_th*0.25*dt*dt*(g/Msq)*(thet-thetm)/(dz*0.5*(thet+thetm)));
                    
                        hz[n] = 0.5 * (hk + hkm) * implicitness * gimp;
                        assert(hz[n] > 0.0); 
                    }
                }
			}
            
            for(k = igz; k < icz - igz; k++) {l = k*icx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j*icx;
                    for(i = igx; i < icx - igx; i++) {n = m + i;
                        hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
                    }
                }
            }
			
			
			break;
		}
		default: ERROR("ndim not in {1,2,3}");
	}
}

/* ========================================================================== */
#ifdef THIRD_ORDER_UPWIND_CORRECTION
#define INTERPOL(l,c,cr) (FIVE_SIXTHS * c + ONE_THIRD * cr - ONE_SIXTH * l)
#endif
#ifdef SECOND_ORDER_CENTRAL_CORRECTION
#define INTERPOL(l,c,cr) (0.5 * (c + cr))
#endif
#ifdef FIRST_ORDER_UPWIND_CORRECTION
#define INTERPOL(l,c,cr) (c)
#endif


static void flux_correction_due_to_pressure_gradients(
                                                      ConsVars* flux[3],
                                                      VectorField* buoy,
                                                      const ElemSpaceDiscr* elem,
                                                      ConsVars* Sol,
                                                      ConsVars* Sol0,
                                                      const MPV* mpv, 
                                                      double* hplus[3],
                                                      double* hgrav,
                                                      double* dp2,
                                                      const double t,
                                                      const double dt,
                                                      const double implicitness) {
    
    extern User_Data ud;
    
    const int ndim = elem->ndim;
    
    int nsp;
    
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
            
            const double dto2dx = implicitness * 0.5 * dt / elem->dx;
            const double dto2dy = implicitness * 0.5 * dt / elem->dy;
            
            ConsVars* f = flux[0];
            ConsVars* g = flux[1];
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            
            double oorhoi, ui, vi, wi, Yi, Zi, Hi, oorhoim, uim, vim, wim, Yim, Zim, Him;
            double Xi[NSPEC], Xim[NSPEC];
            double oorhoj, uj, vj, wj, Yj, Zj, Hj, oorhojm, ujm, vjm, wjm, Yjm, Zjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
            double us, vs, ws, Hs, Ys, Zs;
            double Xs[NSPEC];
            double tmpx, tmpy, frhoY, grhoY;
            
            double p0_c, p0_s;
            double upwind;
            
            int i, j, mc, me, mw, nc, nn, ns, ic, icm, icmm, icp, jc, jcm, jcmm, jcp;
            
            for(j = igy; j < icy - igy; j++) {
                mc    = j * ifx;
                p0_c = mpv->HydroState->p0[j];
                
                for(i = igx; i < ifx - igx; i++) {
                    nn   = mc + i + ifx;
                    nc   = mc + i;
                    ns   = mc + i - ifx;
                    ic   = j*icx + i;
                    icp  = ic + 1;
                    icm  = ic - 1;
                    icmm = ic - 2;
                    
                    buoy->x[ic] = 0.0;
                    
                    /* compute interface values of u, v, Y, Z, H */
                    oorhoi  = 1.0 / INTERPOL(Sol->rhoY[icp],  Sol->rhoY[ic],  Sol->rhoY[icm]); 
                    ui      = oorhoi * INTERPOL(Sol->rhou[icp], Sol->rhou[ic], Sol->rhou[icm]);
                    vi      = oorhoi * INTERPOL(Sol->rhov[icp], Sol->rhov[ic], Sol->rhov[icm]);
                    wi      = oorhoi * INTERPOL(Sol->rhow[icp], Sol->rhow[ic], Sol->rhow[icm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xi[nsp]      = oorhoi * INTERPOL(Sol->rhoX[nsp][icp], Sol->rhoX[nsp][ic], Sol->rhoX[nsp][icm]);
                    }
                    Yi      = oorhoi * INTERPOL(Sol->rho[icp], Sol->rho[ic], Sol->rho[icm]);
                    Zi      = oorhoi * INTERPOL(Sol->rhoZ[icp], Sol->rhoZ[ic], Sol->rhoZ[icm]);
                    Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
                    
                    oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                    uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                    vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                    wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                    }
                    Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                    Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
                    Him     = oorhoim * (INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]) + p0_c);
                    
                    tmpx = - dto2dx * (  0.75  *   hplusx[nc] * (dp2[ic]     - dp2[icm]    )  
                                       + 0.125 * ( hplusx[nn] * (dp2[ic+icx] - dp2[icm+icx])  
                                                 + hplusx[ns] * (dp2[ic-icx] - dp2[icm-icx])  
                                                  ) 
                                       );
                    
                    frhoY = f->rhoY[nc] + tmpx;
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                    
                    us = upwind * uim + (1.0 - upwind) * ui;
                    vs = upwind * vim + (1.0 - upwind) * vi;  
                    ws = upwind * wim + (1.0 - upwind) * wi;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * Xim[nsp] + (1.0 - upwind) * Xi[nsp];
                    }
                    Ys = upwind * Yim + (1.0 - upwind) * Yi;
                    Zs = upwind * Zim + (1.0 - upwind) * Zi;
                    Hs = upwind * Him + (1.0 - upwind) * Hi;
                    
                    f->rho[nc]  = Ys * tmpx;
                    f->rhou[nc] = us * tmpx + us * tmpx;
                    f->rhov[nc] = vs * tmpx;  
                    f->rhow[nc] = ws * tmpx;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        f->rhoX[nsp][nc] = Xs[nsp] * tmpx; 
                    }
                    f->rhoe[nc] = 0.0 * Hs * tmpx;
                    f->rhoY[nc] = tmpx;
                    f->rhoZ[nc] = Zs * tmpx;
                }
            }  
            
            
            /* fluxes in the y-direction */
            for(i = igx; i < icx - igx; i++) {
                nc = i * ify;
                
                for(j = igy; j < ify - igy; j++) {
                    me = nc + j + ify;
                    mc = nc + j;
                    mw = nc + j - ify;
                    jc   = j * icx + i;
                    jcp  = jc + icx;
                    jcm  = jc - icx;
                    jcmm = jc - 2*icx;
                    
                    p0_c = mpv->HydroState->p0[j];
                    p0_s = mpv->HydroState->p0[j-1];
                    
                    buoy->y[jc] = 0.0;
                    
                    oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                    uj = oorhoj * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                    vj = oorhoj * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                    wj = oorhoj * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                    }
                    Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                    Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
                    Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
                    
                    oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                    ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                    vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                    wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                    }
                    Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                    Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
                    Hjm = oorhojm * (INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]) + p0_s);
                    
                    tmpy = - dto2dy * (  0.75  *   hplusy[mc] * (dp2[jc]   - dp2[jcm]  ) 
                                       + 0.125 * ( hplusy[me] * (dp2[jc+1] - dp2[jcm+1]) 
                                                 + hplusy[mw] * (dp2[jc-1] - dp2[jcm-1]) 
                                                  ) 
                                       );
                    
                    tmpy += - 0.5 * dt * hgrav[mc] * 0.5 * (dp2[jc]+dp2[jcm]); /* implicit gravity contribution */                    
                    grhoY       = g->rhoY[mc] + tmpy;
                    
                    
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(grhoY, 0.01) ); 
#endif
                    
                    us = upwind * ujm + (1.0 - upwind) * uj;
                    vs = upwind * vjm + (1.0 - upwind) * vj;
                    ws = upwind * wjm + (1.0 - upwind) * wj;
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * Xjm[nsp] + (1.0 - upwind) * Xj[nsp];
                    }
                    Ys = upwind * Yjm + (1.0 - upwind) * Yj;
                    Zs = upwind * Zjm + (1.0 - upwind) * Zj;
                    Hs = upwind * Hjm + (1.0 - upwind) * Hj;
                    
                    g->rho[mc]  = Ys * tmpy;
                    g->rhou[mc] = us * tmpy;  
                    g->rhov[mc] = vs * tmpy + vs * tmpy; 
                    g->rhow[mc] = ws * tmpy; 
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        g->rhoX[nsp][mc] = Xs[nsp] * tmpy; 
                    }
                    g->rhoe[mc] = 0.0 * Hs * tmpy;
                    g->rhoY[mc] = tmpy;
                    g->rhoZ[mc] = Zs * tmpy;
                }   
            }
            
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
            
            const double dto2dx = implicitness * 0.5 * dt / elem->dx;
            const double dto2dy = implicitness * 0.5 * dt / elem->dy;
            const double dto2dz = implicitness * 0.5 * dt / elem->dz;
            
            const int dix = 1;
            const int diy = icx; 
            const int diz = icx*icy;            
                        
            const double cstencil  = 1.0/64.0;
            
            ConsVars* fx = flux[0];
            ConsVars* fy = flux[1];
            ConsVars* fz = flux[2];
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            const double* hplusz   = hplus[2];
            
            double oorhoi, ui, vi, wi, Yi, Zi, Hi; 
            double oorhoim, uim, vim, wim, Yim, Zim, Him;
            double Xi[NSPEC], Xim[NSPEC];
            
            double oorhoj, uj, vj, wj, Yj, Zj, Hj; 
            double oorhojm, ujm, vjm, wjm, Yjm, Zjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
            
            double oorhok, uk, vk, wk, Yk, Zk, Hk; 
            double oorhokm, ukm, vkm, wkm, Ykm, Zkm, Hkm;
            double Xk[NSPEC], Xkm[NSPEC];
            
            double us, vs, ws, Hs, Ys, Zs;
            double Xs[NSPEC];
            
            double tmpx, tmpy, tmpz, frhoY;
            
            double p0_c, p0_s;
            double upwind;
            
            int i, j, k, l, m, n;
            
            int ic, icm, icmm, icp;
            int jc, jcm, jcmm, jcp;
            int kc, kcm, kcmm, kcp;
            
            for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                    p0_c = mpv->HydroState->p0[j];
                    for(i = igx; i < ifx - igx; i++) {n = m + i;
                        ic   = k*diz + j*diy + i*dix;
                        icp  = ic + dix; /* these appear unused as they only are accessed in macro INTERPOL() */
                        icm  = ic - dix;
                        icmm = ic - 2*dix; /* these appear unused as they only are accessed in macro INTERPOL() */
                        
                        /* compute interface values of u, v, Y, Z, H */
                        
                        oorhoi  = 1.0 / INTERPOL(Sol->rhoY[icp],  Sol->rhoY[ic],  Sol->rhoY[icm]); 
                        ui      = oorhoi * INTERPOL(Sol->rhou[icp], Sol->rhou[ic], Sol->rhou[icm]);
                        vi      = oorhoi * INTERPOL(Sol->rhov[icp], Sol->rhov[ic], Sol->rhov[icm]);
                        wi      = oorhoi * INTERPOL(Sol->rhow[icp], Sol->rhow[ic], Sol->rhow[icm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xi[nsp]      = oorhoi * INTERPOL(Sol->rhoX[nsp][icp], Sol->rhoX[nsp][ic], Sol->rhoX[nsp][icm]);
                        }
                        Yi      = oorhoi * INTERPOL(Sol->rho[icp], Sol->rho[ic], Sol->rho[icm]);
                        Zi      = oorhoi * INTERPOL(Sol->rhoZ[icp], Sol->rhoZ[ic], Sol->rhoZ[icm]);
                        Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
                        
                        oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                        uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                        vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                        wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                        }
                        Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                        Zim     = oorhoim * INTERPOL(Sol->rhoZ[icmm], Sol->rhoZ[icm], Sol->rhoZ[ic]);
                        Him     = oorhoim * (INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]) + p0_c);
                        
                        tmpx = - dto2dx * hplusx[n] * cstencil *
                        (  36.0 *  (dp2[ic] - dp2[icm])  
                         +  6.0 * (  (dp2[ic+diy] - dp2[icm+diy]) 
                                   + (dp2[ic-diy] - dp2[icm-diy])  
                                   + (dp2[ic+diz] - dp2[icm+diz])  
                                   + (dp2[ic-diz] - dp2[icm-diz])
                                   )
                         +  1.0 * (  (dp2[ic+diy+diz] - dp2[icm+diy+diz]) 
                                   + (dp2[ic-diy+diz] - dp2[icm-diy+diz])  
                                   + (dp2[ic+diy-diz] - dp2[icm+diy-diz])  
                                   + (dp2[ic-diy-diz] - dp2[icm-diy-diz])
                                   )
                         );
                        
                        frhoY = fx->rhoY[n] + tmpx;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * uim + (1.0 - upwind) * ui;
                        vs = upwind * vim + (1.0 - upwind) * vi;  
                        ws = upwind * wim + (1.0 - upwind) * wi;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xim[nsp] + (1.0 - upwind) * Xi[nsp];
                        }
                        Ys = upwind * Yim + (1.0 - upwind) * Yi;
                        Zs = upwind * Zim + (1.0 - upwind) * Zi;
                        Hs = upwind * Him + (1.0 - upwind) * Hi;
                        
                        fx->rho[n]  = Ys * tmpx;
                        fx->rhou[n] = us * tmpx + us * tmpx;
                        fx->rhov[n] = vs * tmpx;  
                        fx->rhow[n] = ws * tmpx;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fx->rhoX[nsp][n] = Xs[nsp] * tmpx; 
                        }
                        fx->rhoe[n] = 0.0 * Hs * tmpx;
                        fx->rhoY[n] = tmpx;
                        fx->rhoZ[n] = Zs * tmpx;
                    }
                } 
            }
            
            
            /* fluxes in the y-direction */
            for(i = igx; i < icx - igx; i++) {l = i * ify*icz;
                for(k = igz; k < icz - igz; k++) {m = l + k * ify;
                    for(j = igy; j < ify - igy; j++) {n = m + j;
                        jc   = k*diz + j*diy + i*dix;
                        jcp  = jc + diy;  
                        jcm  = jc - diy;
                        jcmm = jc - 2*diy; 
                        
                        p0_c = mpv->HydroState->p0[j];
                        p0_s = mpv->HydroState->p0[j-1];
                        
                        buoy->y[jc] = 0.0;
                        
                        oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                        uj = oorhoj  * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                        vj = oorhoj  * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                        wj = oorhoj  * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                        }
                        Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                        Zj = oorhoj * INTERPOL(Sol->rhoZ[jcp], Sol->rhoZ[jc], Sol->rhoZ[jcm]);
                        Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
                        
                        oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                        ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                        vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                        wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                        }
                        Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                        Zjm = oorhojm * INTERPOL(Sol->rhoZ[jcmm], Sol->rhoZ[jcm], Sol->rhoZ[jc]);
                        Hjm = oorhojm * (INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]) + p0_s);
                        
                        tmpy = - dto2dy * hplusy[n] * cstencil *
                        (  36.0 *    (dp2[jc] - dp2[jcm])  
                         +  6.0 * (  (dp2[jc+diz] - dp2[jcm+diz]) 
                                   + (dp2[jc-diz] - dp2[jcm-diz])  
                                   + (dp2[jc+dix] - dp2[jcm+dix])  
                                   + (dp2[jc-dix] - dp2[jcm-dix])
                                   )
                         +  1.0 * (  (dp2[jc+diz+dix] - dp2[jcm+diz+dix]) 
                                   + (dp2[jc-diz+dix] - dp2[jcm-diz+dix])  
                                   + (dp2[jc+diz-dix] - dp2[jcm+diz-dix])  
                                   + (dp2[jc-diz-dix] - dp2[jcm-diz-dix])
                                   )
                         );

                        frhoY       = fy->rhoY[n] + tmpy;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * ujm + (1.0 - upwind) * uj;
                        vs = upwind * vjm + (1.0 - upwind) * vj;
                        ws = upwind * wjm + (1.0 - upwind) * wj;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xjm[nsp] + (1.0 - upwind) * Xj[nsp];
                        }
                        Ys = upwind * Yjm + (1.0 - upwind) * Yj;
                        Zs = upwind * Zjm + (1.0 - upwind) * Zj;
                        Hs = upwind * Hjm + (1.0 - upwind) * Hj;
                        
                        fy->rho[n]  = Ys * tmpy;
                        fy->rhou[n] = us * tmpy;  
                        fy->rhov[n] = vs * tmpy + vs * tmpy; 
                        fy->rhow[n] = ws * tmpy;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fy->rhoX[nsp][n] = Xs[nsp] * tmpy; 
                        }
                        fy->rhoe[n] = 0.0 * Hs * tmpy;
                        fy->rhoY[n] = tmpy;
                        fy->rhoZ[n] = Zs * tmpy;
                    }   
                }
            }
            
            /* fluxes in the z-direction */
            for(j = igy; j < icy - igy; j++) {l = j * ifz*icx;
                for(i = igx; i < icx - igx; i++) {m = l + i * ifz;
                    for(k = igz; k < ifz - igz; k++) {n = m + k;
                        kc   = k*diz + j*diy + i*dix;
                        kcp  = kc + diz;   
                        kcm  = kc - diz;
                        kcmm = kc - 2*diz;  
                        
                        p0_c = mpv->HydroState->p0[j];
                        
                        oorhok = 1.0 /INTERPOL(Sol->rhoY[kcp], Sol->rhoY[kc], Sol->rhoY[kcm]);
                        uk = oorhok * INTERPOL(Sol->rhou[kcp], Sol->rhou[kc], Sol->rhou[kcm]);
                        vk = oorhok * INTERPOL(Sol->rhov[kcp], Sol->rhov[kc], Sol->rhov[kcm]);
                        wk = oorhok * INTERPOL(Sol->rhow[kcp], Sol->rhow[kc], Sol->rhow[kcm]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xk[nsp]      = oorhok * INTERPOL(Sol->rhoX[nsp][kcp], Sol->rhoX[nsp][kc], Sol->rhoX[nsp][kcm]);
                        }
                        Yk = oorhok * INTERPOL(Sol->rho[kcp], Sol->rho[kc], Sol->rho[kcm]);
                        Zk = oorhok * INTERPOL(Sol->rhoZ[kcp], Sol->rhoZ[kc], Sol->rhoZ[kcm]);
                        Hk = oorhok * (INTERPOL(Sol->rhoe[kcp], Sol->rhoe[kc], Sol->rhoe[kcm]) + p0_c);
                        
                        oorhokm = 1.0 / INTERPOL(Sol->rhoY[kcmm], Sol->rhoY[kcm], Sol->rhoY[kc]);
                        ukm = oorhokm * INTERPOL(Sol->rhou[kcmm], Sol->rhou[kcm], Sol->rhou[kc]);
                        vkm = oorhokm * INTERPOL(Sol->rhov[kcmm], Sol->rhov[kcm], Sol->rhov[kc]);
                        wkm = oorhokm * INTERPOL(Sol->rhow[kcmm], Sol->rhow[kcm], Sol->rhow[kc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xkm[nsp]      = oorhokm * INTERPOL(Sol->rhoX[nsp][kcmm], Sol->rhoX[nsp][kcm], Sol->rhoX[nsp][kc]);
                        }
                        Ykm = oorhokm * INTERPOL(Sol->rho[kcmm], Sol->rho[kcm], Sol->rho[kc]);
                        Zkm = oorhokm * INTERPOL(Sol->rhoZ[kcmm], Sol->rhoZ[kcm], Sol->rhoZ[kc]);
                        Hkm = oorhokm * (INTERPOL(Sol->rhoe[kcmm], Sol->rhoe[kcm], Sol->rhoe[kc]) + p0_c);
                        
                        tmpz = - dto2dz * hplusz[n] * cstencil *
                        (  36.0 *    (dp2[kc] - dp2[kcm])  
                         +  6.0 * (  (dp2[kc+diz] - dp2[kcm+diz]) 
                                   + (dp2[kc-diz] - dp2[kcm-diz])  
                                   + (dp2[kc+dix] - dp2[kcm+dix])  
                                   + (dp2[kc-dix] - dp2[kcm-dix])
                                   )
                         +  1.0 * (  (dp2[kc+diz+dix] - dp2[kcm+diz+dix]) 
                                   + (dp2[kc-diz+dix] - dp2[kcm-diz+dix])  
                                   + (dp2[kc+diz-dix] - dp2[kcm+diz-dix])  
                                   + (dp2[kc-diz-dix] - dp2[kcm-diz-dix])
                                   )
                         );

                        frhoY       = fz->rhoY[n] + tmpz;
                        
#ifdef NO_UPWIND_PROJ1
                        upwind = 0.5; 
#else
                        upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                        
                        us = upwind * ukm + (1.0 - upwind) * uk;
                        vs = upwind * vkm + (1.0 - upwind) * vk;
                        ws = upwind * wkm + (1.0 - upwind) * wk;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xs[nsp]      = upwind * Xkm[nsp] + (1.0 - upwind) * Xk[nsp];
                        }
                        Ys = upwind * Ykm + (1.0 - upwind) * Yk;
                        Zs = upwind * Zkm + (1.0 - upwind) * Zk;
                        Hs = upwind * Hkm + (1.0 - upwind) * Hk;
                        
                        fz->rho[n]  = Ys * tmpz;
                        fz->rhou[n] = us * tmpz;  
                        fz->rhov[n] = vs * tmpz; 
                        fz->rhow[n] = ws * tmpz + ws * tmpz;  
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            fz->rhoX[nsp][n] = Xs[nsp] * tmpz; 
                        }
                        fz->rhoe[n] = 0.0 * Hs * tmpz;
                        fz->rhoY[n] = tmpz;
                        fz->rhoZ[n] = Zs * tmpz;
                    }   
                }
            }			
            break;
        }
        default: ERROR("ndim not in {1, 2, 3}");
    }
}

/* ========================================================================== */

static void rhs_fix_for_open_boundaries(
										double* rhs, 
										const ElemSpaceDiscr* elem, 
										ConsVars* Sol_new, 
										ConsVars* Sol_old, 
										ConsVars* flux[3],
										double dt,
										MPV* mpv) {
	
	extern User_Data ud;
	
	const States* HydroState = mpv->HydroState;
	const double implicitness = 0.5*ud.implicitness;
	
	double** hplus       = mpv->Level[0]->wplus;
	
	int ndim = elem->ndim;
	
	if (elem->left == OPEN) {
		
		assert(elem->left == elem->right); /* if you get thrown out here: one-sided open domain not yet impltd */
		
		switch(ndim) {
			case 1: {
				ERROR("function not available");
				break;
			}
				
			case 2: {
				
				const double factor = 1.0 / (implicitness * elem->dx * dt);
				
				const int igx = elem->igx;
				const int icx = elem->icx;
				const int ifx = elem->ifx;
				
				const int igy = elem->igy;
				const int icy = elem->icy;
				
				double* hplusx  = hplus[0];
				
				double rhs_sum, coeff_sum_right, coeff_sum_left, du;
				int j, mc, nc_left, nc_right, mf, nf_left, nf_right, count;
				
				count = 0;
				rhs_sum  = 0.0;
				coeff_sum_left  = 0.0;
				coeff_sum_right = 0.0;
				
				for(j = igy; j < icy - igy; j++) {
					mf = j * ifx;
					nf_left  = mf + igx;
					nf_right = mf + ifx - igx - 1;
					
					rhs_sum         += (flux[0]->rhoe[nf_right] - flux[0]->rhoe[nf_left]);
					coeff_sum_left  += HydroState->rho0[j] * hplusx[nf_left]; 
					coeff_sum_right += HydroState->rho0[j] * hplusx[nf_right];
					
					count++; 
				} 
				
				rhs_sum         /= count;
				coeff_sum_left  /= count;
				coeff_sum_right /= count;
				du = - rhs_sum / (coeff_sum_left + coeff_sum_right);
				
				for(j = igy; j < icy - igy; j++) {
					mc = j * icx;
					mf = j * ifx;
					nc_left  = mc + igx;
					nc_right = mc + icx - igx - 1;
					nf_left  = mf + igx;
					nf_right = mf + ifx - igx - 1;
					
					rhs[nc_left]  -= factor * du * HydroState->rho0[j] * hplusx[nf_left];
					rhs[nc_right] += factor * du * HydroState->rho0[j] * hplusx[nf_right];
					
				}
				
				mpv->du = du;
				
				break;
			}
			case 3: {
				ERROR("ndim = 3 - option not implemented (kinetic_energy_change_explicit() )");
			}
			default: ERROR("ndim not in {1, 2, 3}");
		}
	}
}

/* ========================================================================== */

static void flux_correction_due_to_pressure_values(
												   ConsVars* flux[3],
												   VectorField* buoy,
												   const ElemSpaceDiscr* elem,
												   ConsVars* Sol,
												   double* dp2, 
												   const double dt) {
	
	extern User_Data ud;
    extern Thermodynamic th;
	
    const double Gammainv = th.Gammainv;
    
	const int ndim = elem->ndim;
	    
	switch(ndim) {
		case 1: {
			ERROR("function not available");
			break;
		}
			
		case 2: {
            double coeff;
			int i, j, m, n, ic, icm, jc, jcm;
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int ifx = elem->ifx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int ify = elem->ify;
			
			ConsVars* f = flux[0];
			ConsVars* g = flux[1];
			
            /* flux_post_correction    change from   flux ... +=  corr    to    flux ... = corr */
			
			/* Momentum flux corrections */
			for(j = igy; j < icy - igy; j++) {
				
				m = j * ifx;
				
				for(i = igx; i < ifx - igx; i++) {
					
					n   = m + i;
					ic  = n - j;
					icm = ic - 1;
					
#ifdef THERMCON
                    coeff = Gammainv * 0.5 * (Sol->rhoY[ic]*Sol->rhoY[ic]/Sol->rho[ic] + Sol->rhoY[icm]*Sol->rhoY[icm]/Sol->rho[icm]);
#else
                    coeff = 1.0;
#endif
                    f->rhou[n] += ud.p_extrapol*0.5*coeff*(
                                              (ud.latw[0]*dp2[ic+icx]  + ud.latw[1]*dp2[ic]  + ud.latw[2]*dp2[ic-icx]) 
                                            + (ud.latw[0]*dp2[icm+icx] + ud.latw[1]*dp2[icm] + ud.latw[2]*dp2[icm-icx]) 
                                               );
				}
			}  
			
			/* fluxes in the y-direction */
			for(i = igx; i < icx - igx; i++) {
				
				n = i * ify;
				
				for(j = igy; j < ify - igy; j++) {
					
					m   = n + j;
					jc  = j * icx + i;
					jcm = jc - icx;
					
#ifdef THERMCON
                    coeff = Gammainv * 0.5 * (Sol->rhoY[jc]*Sol->rhoY[jc]/Sol->rho[jc] + Sol->rhoY[jcm]*Sol->rhoY[jcm]/Sol->rho[jcm]);
#else
                    coeff = 1.0;
#endif
					g->rhov[m] += ud.p_extrapol*0.5*coeff*(  
                                                     (ud.latw[0]*dp2[jc+1]  + ud.latw[1]*dp2[jc]  + ud.latw[2]*dp2[jc-1])  
                                                   + (ud.latw[0]*dp2[jcm+1] + ud.latw[1]*dp2[jcm] + ud.latw[2]*dp2[jcm-1])
                                                  );
				}   
			}
			
			break;
		}
		case 3: {
			ERROR("function not available");
			break;
		}
		default: ERROR("ndim not in {1, 2,3}");
    }
#if 0
    printf("leaving flux_correction_due_to_pressure_values()\n");
#endif
}

/* ========================================================================== */

static void flux_fix_for_open_boundaries(
										 ConsVars* flux[3],
										 const ElemSpaceDiscr* elem, 
										 MPV* mpv) {
	
	const States* HydroState = mpv->HydroState;    
	double** hplus           = mpv->Level[0]->wplus;
	double du                = mpv->du;
	double Sdu               = mpv->Sdu;
	
	int ndim = elem->ndim;
	
	if (elem->left == OPEN) {
		
		assert(elem->left == elem->right); /* if you get thrown out here: one-sided open domain not yet impltd */
		
		switch(ndim) {
			case 1: {
				ERROR("function not available");
				break;
			}
				
			case 2: {
				
				const int igx = elem->igx;
				const int ifx = elem->ifx;
				
				const int igy = elem->igy;
				const int icy = elem->icy;
				
				double* hplusx  = hplus[0];
				
				int j, mf, nf_left, nf_right;
				
				for(j = igy; j < icy - igy; j++) {
					mf = j * ifx;
					nf_left  = mf + igx;
					nf_right = mf + ifx - igx - 1;
					
					flux[0]->rho[nf_left]   -= Sdu * HydroState->rho0[j];
					flux[0]->rho[nf_right]  += Sdu * HydroState->rho0[j];
					flux[0]->rhoe[nf_left]  -= du  * HydroState->rho0[j] * hplusx[nf_left];
					flux[0]->rhoe[nf_right] += du  * HydroState->rho0[j] * hplusx[nf_right];
				} 
				
				break;
			}
			case 3: {
				ERROR("ndim = 3 - option not implemented (kinetic_energy_change_explicit() )");
			}
			default: ERROR("ndim not in {1, 2, 3}");
		}
	}
}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
