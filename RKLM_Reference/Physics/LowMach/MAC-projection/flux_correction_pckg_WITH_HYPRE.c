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

static enum Constraint integral_condition(ConsVars* flux[3],
										  double* rhs, 
										  ConsVars* Sol,
										  const double dt,
										  const ElemSpaceDiscr* elem,
										  MPV* mpv);

static void controlled_variable_change_explicit(double* rhs, 
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
                                                      double* hS,
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

void flux_correction(ConsVars* flux[3],
                     VectorField* adv_flux_diff,
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
    double*  hS          = mpv->Level[0]->wgrav;
	
	double* rhs          = mpv->Level[0]->rhs;
	double* dp2		     = mpv->Level[0]->p;
        
    double rhsmax;
    
	int n;
        
    printf("\n\n====================================================");
    printf("\nFirst Projection");
    printf("\n====================================================\n");
	    
	operator_coefficients(hplus, hcenter, hS, elem, Sol, Sol0, mpv, dt);

    if (ud.time_integrator == SI_MIDPT) {
        rhsmax = controlled_variable_flux_divergence(rhs, (const ConsVars**)flux, dt, elem);
        assert(integral_condition(flux, rhs, Sol, dt, elem, mpv) != VIOLATED); 
        if (ud.is_compressible) {
            for (int nc=0; nc<elem->nc; nc++) {
                rhs[nc] += hcenter[nc]*mpv->p2_cells[nc];
            }
        }
    } else {
        controlled_variable_change_explicit(rhs, elem, Sol, Sol0, dt, mpv);        
        assert(integral_condition(flux, rhs, Sol, dt, elem, mpv) != VIOLATED); 
    }
    
    rhs_fix_for_open_boundaries(rhs, elem, Sol, Sol0, flux, dt, mpv);
    
#if 0
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
    
    sprintf(fn2, "%s/Tests/frhoY_y_pre.hdf", ud.file_name);
    sprintf(fieldname2, "frhoY_y_pre.hdf");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->ify,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             flux[1]->rhoY,
             fn2,
             fieldname2);

    /*
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
    variable_coefficient_poisson_cells(dp2, rhs, (const double **)hplus, hcenter, Sol, elem, node, implicitness);
    set_ghostcells_p2(dp2, (const double **)hplus, elem, elem->igx);

    /* Note: flux will contain only the flux-correction after this routine; 
     it is thus overwritten under the SI_MIDPT time integration sequence */
    flux_correction_due_to_pressure_gradients(flux, buoy, elem, Sol, Sol0, mpv, hplus, hS, dp2, t, dt, implicitness);
    
#if 0
    sprintf(fn2, "%s/Tests/frhoY_y_post.hdf", ud.file_name);
    sprintf(fieldname2, "frhoY_y_post.hdf");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->ify,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             flux[1]->rhoY,
             fn2,
             fieldname2);
#endif
    
    /* under the SI_MIDPT time integrator sequence, after this step flux will contain:
       (flux_correction + recomputed flux after first advection sequence)
        - flux accumulated over first advection sequence
       The sum in the bracket is a divergence-controlled flux, while the 
       subtracted flux from the first sequence will make up for the ``bad guess''
       of the flux divergence in the explicity advection cycle.
     */
    if (ud.time_integrator == SI_MIDPT) {
        for (int ii=0; ii<elem->nfx; ii++) flux[0]->rhoY[ii] += adv_flux_diff->x[ii];
        if (elem->ndim > 1) for (int ii=0; ii<elem->nfy; ii++) flux[1]->rhoY[ii] += adv_flux_diff->y[ii];
        if (elem->ndim > 2) for (int ii=0; ii<elem->nfz; ii++) flux[2]->rhoY[ii] += adv_flux_diff->z[ii];
    }
    if (ud.p_flux_correction) {
        flux_correction_due_to_pressure_values(flux, buoy, elem, Sol, dp2, dt); 
    }

#if 0
    sprintf(fn2, "%s/Tests/frhoY_y_post2.hdf", ud.file_name);
    sprintf(fieldname2, "frhoY_y_post2.hdf");
    
    WriteHDF(prhsfile,
             mpv->Level[0]->elem->ify,
             mpv->Level[0]->elem->icx,
             mpv->Level[0]->elem->icz,
             mpv->Level[0]->elem->ndim,
             flux[1]->rhoY,
             fn2,
             fieldname2);
#endif

    
    flux_fix_for_open_boundaries(flux, elem, mpv);  

    memcpy(mpv->dp2_cells, dp2, elem->nc*sizeof(double));

    /* store results in mpv-fields */
    if (ud.time_integrator == SI_MIDPT) {
        for(n=0; n<elem->nc; n++) {
            /* goal is to actually compute full deviation from background pressure here to begin with */
            mpv->dp2_cells[n] = mpv->p2_cells[n] = dp2[n];
        }
    } else {
        for(n=0; n<elem->nc; n++) {
            mpv->dp2_cells[n] = Sol->rhoZ[PRES][n] + (1.0-DP2_ELL_FACTOR)*dp2[n] - mpv->p2_cells[n];
            mpv->p2_cells[n]  = mpv->p2_cells[n] + mpv->dp2_cells[n];
        }        
    }
    
	set_ghostcells_p2(mpv->p2_cells, (const double **)hplus, elem, elem->igx);
}



/* ========================================================================== */

double controlled_variable_flux_divergence(double* rhs, 
                                           const ConsVars* flux[3],
                                           const double dt, 
                                           const ElemSpaceDiscr* elem)
{
    /* right hand side of pressure equation via advective flux divergence */
    memset(rhs, 0.0, elem->nc*sizeof(double));
    
    const int igx = elem->igx;
    const int icx = elem->icx;
    const int ifx = elem->ifx;
    const int igy = elem->igy;
    const int icy = elem->icy;
    const int ify = elem->ify;
    
    const double factor = 2.0 / dt;
    const double dx = elem->dx;
    const double dy = elem->dy;
    
    double rhsmax = 0.0;
    
    assert(elem->ndim == 2);
    
    for(int i=0; i<elem->nc; i++) rhs[i] = 0.0;
    
    for(int j = igy; j < icy - igy; j++) {
        int mc  = j * icx;
        int mfx = j * ifx; 
        int mfy = j;
        for(int i = igx; i < icx - igx; i++) {
            int nc  = mc  + i;
            int nfx = mfx + i;
            int nfy = mfy + i*ify;
            rhs[nc] = factor * ((flux[0]->rhoY[nfx+1] - flux[0]->rhoY[nfx])/dx + (flux[1]->rhoY[nfy+1] - flux[1]->rhoY[nfy])/dy);
            /* rhs[nc]+= hcenter[nc]*mpv->p2_cells[nc]; */
            rhsmax  = MAX_own(rhsmax, fabs(rhs[nc]));
        }
    }
    
#if 0
    FILE *prhsfile = NULL;
    char fn[100], fieldname[90];
    sprintf(fn, "rhs_cells.hdf");
    sprintf(fieldname, "rhs-cells");    
    WriteHDF(prhsfile, elem->icx, elem->icy, elem->icz, elem->ndim, rhs, fn, fieldname);
#endif
    
    return rhsmax;
}


/* ========================================================================== */
 
static void controlled_variable_change_explicit(double* rhs, 
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
                           double* hS,
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
    
	const double implicitness = ud.implicitness;
    const double ccw = (ud.time_integrator == SI_MIDPT ? 4.0 : 2.0); 
    /* when p2 is perturbation pressure:
     const double ccenter = - ccw * (ud.compressibility*ud.Msq)*th.gamminv/(mpv->dt*mpv->dt); 
     const double cexp    = 1.0-th.gamm;
     */
    /* when p2 is Exner pressure: */
    const double ccenter = - ccw * (ud.compressibility*ud.Msq)*th.gm1inv/(mpv->dt*mpv->dt); 
    const double cexp    = 2.0-th.gamm;
	
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
            const double dy = elem->dy;
			double* hx = hplus[0];
			double* hy = hplus[1];
			double* hc = wcenter;
			
			double hi, him, hj, hjm, g, gimp, Msq;
			
			int i, j, m, n, ic, icm, jc, jcm;
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[0];
            
			for(j = igy-1; j < icy - igy+1; j++) {
                m = j * ifx;
				
                for(i = igx-1; i < ifx - igx+1; i++) {
                    n     = m + i;
                    ic    = n - j;
					icm   = ic - 1; 
					
#ifdef TIME_AVERAGED_COEFFS_PROJ1
                    hi    = 0.5 * ( Sol->rhoY[ic] * Sol->rhoY[ic] / Sol->rho[ic] + Sol0->rhoY[ic] * Sol0->rhoY[ic] / Sol0->rho[ic] ) * Gammainv;   
                    him   = 0.5 * ( Sol->rhoY[icm] * Sol->rhoY[icm] / Sol->rho[icm] + Sol0->rhoY[icm] * Sol0->rhoY[icm] / Sol0->rho[icm] ) * Gammainv;
#else
                    hi    = Sol->rhoY[ic] * Sol->rhoY[ic] / Sol->rho[ic]    * Gammainv;   
                    him   = Sol->rhoY[icm] * Sol->rhoY[icm] / Sol->rho[icm] * Gammainv;
#endif
                    
					hx[n] = 0.5 * (hi + him) * implicitness;
                    
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
				
#ifdef TIME_AVERAGED_COEFFS_PROJ1
                    hj    = 0.5 * ( Sol->rhoY[jc] * Sol->rhoY[jc] / Sol->rho[jc] + Sol0->rhoY[jc] * Sol0->rhoY[jc] / Sol0->rho[jc]) * Gammainv;
                    hjm   = 0.5 * ( Sol->rhoY[jcm] * Sol->rhoY[jcm] / Sol->rho[jcm] + Sol0->rhoY[jcm] * Sol0->rhoY[jcm] / Sol0->rho[jcm]) * Gammainv;
#else
                    hj    = Sol->rhoY[jc] * Sol->rhoY[jc] / Sol->rho[jc] * Gammainv;
                    hjm   = Sol->rhoY[jcm] * Sol->rhoY[jcm] / Sol->rho[jcm] * Gammainv;
#endif
                    
                    double S     = mpv->HydroState->S0[j];
                    double Sm    = mpv->HydroState->S0[j-1];
                    double Y     = 0.5 * (Sol->rhoY[jc]  / Sol->rho[jc]  + Sol0->rhoY[jc]  / Sol0->rho[jc]);
                    double Ym    = 0.5 * (Sol->rhoY[jcm] / Sol->rho[jcm] + Sol0->rhoY[jcm] / Sol0->rho[jcm]);
                    double Nsq   = - (g/Msq) * 0.5*(Y+Ym) * (S-Sm)/dy;
                    /* double Nsqsc = impl_grav_th*0.5*dt*dt*Nsq; */
                    double Nsqsc = impl_grav_th*0.25*dt*dt*Nsq;
                     
                    gimp  = 1.0 / (1.0 + Nsqsc);
                    
					hy[m] = 0.5 * (hj + hjm) * implicitness * gimp;
                    hS[m] = Nsqsc * gimp * Gammainv / (g/Msq) / Y;
                                        
					assert(hy[m] > 0.0);
				}
			}
            
			for(j = igy; j < icy - igy; j++) {m = j * icx;
				for(i = igx; i < icx - igx; i++) {n = m + i;
#ifdef TIME_AVERAGED_COEFFS_PROJ1
					hc[n] = ccenter * pow(0.5*(Sol->rhoY[n]+Sol0->rhoY[n]),cexp);
#else
                    hc[n] = ccenter * pow(Sol0->rhoY[n],cexp);
#endif
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
            const double dy = elem->dy;
			double* hx = hplus[0];
			double* hy = hplus[1];
			double* hz = hplus[2];
			double* hc = wcenter;
			
			double hi, him, hj, hjm, hk, hkm, g, gimp, Msq;
			
			int i, j, k, l, m, n, ic, icm, jc, jcm, kc, kcm;
            
            Msq = ud.Msq;
            g   = ud.gravity_strength[0];
            assert(g==0.0); /* implicit gravity only for y-direction */
            
			for(k = igz; k < icz - igz; k++) {l = k * ifx*icy;
                for(j = igy; j < icy - igy; j++) {m = l + j * ifx;
                    for(i = igx; i < ifx - igx; i++) {n = m + i;
                        ic  = k*icx*icy + j*icx + i;
                        icm = ic - 1; 

                        hi    = 0.5 * (Sol->rhoY[ic] *Sol->rhoY[ic] /Sol->rho[ic]  + Sol0->rhoY[ic] *Sol0->rhoY[ic] /Sol0->rho[ic] ) * Gammainv;   
                        him   = 0.5 * (Sol->rhoY[icm]*Sol->rhoY[icm]/Sol->rho[icm] + Sol0->rhoY[icm]*Sol0->rhoY[icm]/Sol0->rho[icm]) * Gammainv;
					
                        /* optional with new time level data only as in 2D part of routine ... */

                        hx[n] = 0.5 * (hi + him) * implicitness;
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

                        hj       = 0.5 * (Sol->rhoY[jc] *Sol->rhoY[jc] /Sol->rho[jc] + Sol0->rhoY[jc] *Sol0->rhoY[jc] /Sol0->rho[jc] ) * Gammainv;
                        hjm      = 0.5 * (Sol->rhoY[jcm]*Sol->rhoY[jcm]/Sol->rho[jcm]+ Sol0->rhoY[jcm]*Sol0->rhoY[jcm]/Sol0->rho[jcm]) * Gammainv;

                        /* optional with new time level data only as in 2D part of routine ... */

                        double S   = mpv->HydroState->S0[j];
                        double Sm  = mpv->HydroState->S0[j-1];
                        double Y   = 0.5 * (Sol->rhoY[jc]  / Sol->rho[jc]  + Sol0->rhoY[jc]  / Sol0->rho[jc]);
                        double Ym  = 0.5 * (Sol->rhoY[jcm] / Sol->rho[jcm] + Sol0->rhoY[jcm] / Sol0->rho[jcm]);
                        double Nsq = - (g/Msq) * 0.5*(Y+Ym) * (S-Sm)/dy;
                        
                        gimp  = 1.0 / (1.0 + impl_grav_th*0.5*dt*dt*Nsq);
                    
                        hy[n] = 0.5 * (hj + hjm) * implicitness * gimp;
                        
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

                        hk       = 0.5 * (Sol->rhoY[kc] *Sol->rhoY[kc] /Sol->rho[kc] + Sol0->rhoY[kc] *Sol0->rhoY[kc] /Sol0->rho[kc] ) * Gammainv;
                        hkm      = 0.5 * (Sol->rhoY[kcm]*Sol->rhoY[kcm]/Sol->rho[kcm]+ Sol0->rhoY[kcm]*Sol0->rhoY[kcm]/Sol0->rho[kcm]) * Gammainv;

                        /* optional with new time level data only as in 2D part of routine ... */
                    
                        hz[n] = 0.5 * (hk + hkm) * implicitness;
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
                                                      double* hS,
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
            
            const double w  = 1.0;
            const double w0 = 1.0-w;
            
            const int igx = elem->igx;
            const int icx = elem->icx;
            const int ifx = elem->ifx;
            const int igy = elem->igy;
            const int icy = elem->icy;
            const int ify = elem->ify;
            
            const double dto2dx = implicitness * 0.5 * dt / elem->dx;
            const double dto2dy = implicitness * 0.5 * dt / elem->dy;
            const double oody   = 1.0/elem->dy;
            
            ConsVars* f = flux[0];
            ConsVars* g = flux[1];
            
            const double* hplusx   = hplus[0];
            const double* hplusy   = hplus[1];
            
            double oorhoi, ui, vi, wi, Yi, Hi, oorhoim, uim, vim, wim, Yim, Him;
            double Xi[NSPEC], Xim[NSPEC];
            double oorhoj, uj, vj, wj, Yj, Hj, oorhojm, ujm, vjm, wjm, Yjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
            double us, vs, ws, Hs, Ys;
            double Xs[NSPEC];

            double oorhoi0, ui0, vi0, wi0, Yi0, Hi0, oorhoim0, uim0, vim0, wim0, Yim0, Him0;
            double Xi0[NSPEC], Xim0[NSPEC];
            double oorhoj0, uj0, vj0, wj0, Yj0, Hj0, oorhojm0, ujm0, vjm0, wjm0, Yjm0, Hjm0;
            double Xj0[NSPEC], Xjm0[NSPEC];
            double us0, vs0, ws0, Hs0, Ys0;
            double Xs0[NSPEC];

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
                    
                    /* TODO: Check what happens if you interpolate (rhoX / rhoY) directly instead */
                    /* compute interface values of u, v, Y, Z, H */
                    oorhoi  = 1.0 / INTERPOL(Sol->rhoY[icp],  Sol->rhoY[ic],  Sol->rhoY[icm]); 
                    ui      = oorhoi * INTERPOL(Sol->rhou[icp], Sol->rhou[ic], Sol->rhou[icm]);
                    vi      = oorhoi * INTERPOL(Sol->rhov[icp], Sol->rhov[ic], Sol->rhov[icm]);
                    wi      = oorhoi * INTERPOL(Sol->rhow[icp], Sol->rhow[ic], Sol->rhow[icm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xi[nsp]      = oorhoi * INTERPOL(Sol->rhoX[nsp][icp], Sol->rhoX[nsp][ic], Sol->rhoX[nsp][icm]);
                    }
                    Yi      = oorhoi * INTERPOL(Sol->rho[icp], Sol->rho[ic], Sol->rho[icm]);
                    Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
                    
                    oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                    uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                    vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                    wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                    }
                    Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
                    Him     = oorhoim * (INTERPOL(Sol->rhoe[icmm], Sol->rhoe[icm], Sol->rhoe[ic]) + p0_c);
          
                    oorhoi0  = 1.0 / INTERPOL(Sol0->rhoY[icp],  Sol0->rhoY[ic],  Sol0->rhoY[icm]); 
                    ui0      = oorhoi0 * INTERPOL(Sol0->rhou[icp], Sol0->rhou[ic], Sol0->rhou[icm]);
                    vi0      = oorhoi0 * INTERPOL(Sol0->rhov[icp], Sol0->rhov[ic], Sol0->rhov[icm]);
                    wi0      = oorhoi0 * INTERPOL(Sol0->rhow[icp], Sol0->rhow[ic], Sol0->rhow[icm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xi0[nsp]      = oorhoi0 * INTERPOL(Sol0->rhoX[nsp][icp], Sol0->rhoX[nsp][ic], Sol0->rhoX[nsp][icm]);
                    }
                    Yi0      = oorhoi0 * INTERPOL(Sol0->rho[icp], Sol0->rho[ic], Sol0->rho[icm]);
                    Hi0      = oorhoi0 * (INTERPOL(Sol0->rhoe[icp], Sol0->rhoe[ic], Sol0->rhoe[icm]) + p0_c);
                    
                    oorhoim0 = 1.0 / INTERPOL(Sol0->rhoY[icmm], Sol0->rhoY[icm], Sol0->rhoY[ic]); 
                    uim0     = oorhoim0 * INTERPOL(Sol0->rhou[icmm], Sol0->rhou[icm], Sol0->rhou[ic]);
                    vim0     = oorhoim0 * INTERPOL(Sol0->rhov[icmm], Sol0->rhov[icm], Sol0->rhov[ic]);
                    wim0     = oorhoim0 * INTERPOL(Sol0->rhow[icmm], Sol0->rhow[icm], Sol0->rhow[ic]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xim0[nsp]      = oorhoim0 * INTERPOL(Sol0->rhoX[nsp][icmm], Sol0->rhoX[nsp][icm], Sol0->rhoX[nsp][ic]);
                    }
                    Yim0     = oorhoim0 * INTERPOL(Sol0->rho[icmm], Sol0->rho[icm], Sol0->rho[ic]);
                    Him0     = oorhoim0 * (INTERPOL(Sol0->rhoe[icmm], Sol0->rhoe[icm], Sol0->rhoe[ic]) + p0_c);

                    tmpx = - dto2dx * (  0.75  *   hplusx[nc] * (dp2[ic]     - dp2[icm]    )  
                                       + 0.125 * ( hplusx[nn] * (dp2[ic+icx] - dp2[icm+icx])  
                                                 + hplusx[ns] * (dp2[ic-icx] - dp2[icm-icx])  
                                                  ) 
                                       );
                    
                    frhoY = f->rhoY[nc] + tmpx;
#ifdef CENTERED_UPDATE_1ST_PROJ
                    buoy->x[ic]  += tmpx; 
                    buoy->x[icm] += tmpx; 
#endif
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(frhoY, 0.01) ); 
#endif
                    
                    us = upwind * (w*uim + w0*uim0) + (1.0 - upwind) * (w*ui + w0*ui0);
                    vs = upwind * (w*vim + w0*vim0) + (1.0 - upwind) * (w*vi + w0*vi0);  
                    ws = upwind * (w*wim + w0*wim0) + (1.0 - upwind) * (w*wi + w0*wi0);  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * (w*Xim[nsp] + w0*Xim0[nsp]) + (1.0 - upwind) * (w*Xi[nsp] + w0*Xi0[nsp]);
                    }
                    Ys = upwind * (w*Yim + w0*Yim0) + (1.0 - upwind) * (w*Yi + w0*Yi0);
                    Hs = upwind * (w*Him + w0*Him0) + (1.0 - upwind) * (w*Hi + w0*Hi0);
                    
                    f->rho[nc]  = Ys * tmpx;
                    f->rhou[nc] = us * tmpx + us * tmpx;
                    f->rhov[nc] = vs * tmpx;  
                    f->rhow[nc] = ws * tmpx;  
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        f->rhoX[nsp][nc] = Xs[nsp] * tmpx; 
                    }
                    f->rhoe[nc] = 0.0 * Hs * tmpx;
                    f->rhoY[nc] = tmpx;
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
                    
                    /* TODO: eliminate buoy-fields from projection step(s); they are everywhere zero anyway it seems. */
                    buoy->y[jc] = 0.0;
                    
                    oorhoj = 1.0 / INTERPOL(Sol->rhoY[jcp], Sol->rhoY[jc], Sol->rhoY[jcm]);
                    uj = oorhoj * INTERPOL(Sol->rhou[jcp], Sol->rhou[jc], Sol->rhou[jcm]);
                    vj = oorhoj * INTERPOL(Sol->rhov[jcp], Sol->rhov[jc], Sol->rhov[jcm]);
                    wj = oorhoj * INTERPOL(Sol->rhow[jcp], Sol->rhow[jc], Sol->rhow[jcm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xj[nsp]      = oorhoj * INTERPOL(Sol->rhoX[nsp][jcp], Sol->rhoX[nsp][jc], Sol->rhoX[nsp][jcm]);
                    }
                    Yj = oorhoj * INTERPOL(Sol->rho[jcp], Sol->rho[jc], Sol->rho[jcm]);
                    Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
                    
                    oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                    ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                    vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                    wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                    }
                    Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
                    Hjm = oorhojm * (INTERPOL(Sol->rhoe[jcmm], Sol->rhoe[jcm], Sol->rhoe[jc]) + p0_s);

                    oorhoj0 = 1.0 / INTERPOL(Sol0->rhoY[jcp], Sol0->rhoY[jc], Sol0->rhoY[jcm]);
                    uj0 = oorhoj0 * INTERPOL(Sol0->rhou[jcp], Sol0->rhou[jc], Sol0->rhou[jcm]);
                    vj0 = oorhoj0 * INTERPOL(Sol0->rhov[jcp], Sol0->rhov[jc], Sol0->rhov[jcm]);
                    wj0 = oorhoj0 * INTERPOL(Sol0->rhow[jcp], Sol0->rhow[jc], Sol0->rhow[jcm]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xj0[nsp]      = oorhoj0 * INTERPOL(Sol0->rhoX[nsp][jcp], Sol0->rhoX[nsp][jc], Sol0->rhoX[nsp][jcm]);
                    }
                    Yj0 = oorhoj0 * INTERPOL(Sol0->rho[jcp], Sol0->rho[jc], Sol0->rho[jcm]);
                    Hj0 = oorhoj0 * (INTERPOL(Sol0->rhoe[jcp], Sol0->rhoe[jc], Sol0->rhoe[jcm]) + p0_c);
                    
                    oorhojm0 = 1.0 / INTERPOL(Sol0->rhoY[jcmm], Sol0->rhoY[jcm], Sol0->rhoY[jc]);
                    ujm0 = oorhojm0 * INTERPOL(Sol0->rhou[jcmm], Sol0->rhou[jcm], Sol0->rhou[jc]);
                    vjm0 = oorhojm0 * INTERPOL(Sol0->rhov[jcmm], Sol0->rhov[jcm], Sol0->rhov[jc]);
                    wjm0 = oorhojm0 * INTERPOL(Sol0->rhow[jcmm], Sol0->rhow[jcm], Sol0->rhow[jc]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xjm0[nsp]      = oorhojm0 * INTERPOL(Sol0->rhoX[nsp][jcmm], Sol0->rhoX[nsp][jcm], Sol0->rhoX[nsp][jc]);
                    }
                    Yjm0 = oorhojm0 * INTERPOL(Sol0->rho[jcmm], Sol0->rho[jcm], Sol0->rho[jc]);
                    Hjm0 = oorhojm0 * (INTERPOL(Sol0->rhoe[jcmm], Sol0->rhoe[jcm], Sol0->rhoe[jc]) + p0_s);

                    tmpy = - dto2dy * (  0.75  *   hplusy[mc] * (dp2[jc]   - dp2[jcm]  ) 
                                       + 0.125 * ( hplusy[me] * (dp2[jc+1] - dp2[jcm+1]) 
                                                 + hplusy[mw] * (dp2[jc-1] - dp2[jcm-1]) 
                                                  ) 
                                       );
                    
                    grhoY        = g->rhoY[mc] + tmpy;
#ifdef CENTERED_UPDATE_1ST_PROJ
                    buoy->y[jc]  += tmpy;
                    buoy->y[jcm] += tmpy;
#endif 
                    
#ifdef NO_UPWIND_PROJ1
                    upwind = 0.5; 
#else
                    upwind = 0.5 * ( 1.0 + SMOOTHSIGN(grhoY, 0.01) ); 
#endif
                    
                    us = upwind * (w*ujm + w0*ujm0) + (1.0 - upwind) * (w*uj + w0*uj0);
                    vs = upwind * (w*vjm + w0*vjm0) + (1.0 - upwind) * (w*vj + w0*vj0);
                    ws = upwind * (w*wjm + w0*wjm0) + (1.0 - upwind) * (w*wj + w0*wj0);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        Xs[nsp]      = upwind * (w*Xjm[nsp] + w0*Xjm0[nsp]) + (1.0 - upwind) * (w*Xj[nsp] + w0*Xj0[nsp]);
                    }
                    Ys = upwind * (w*Yjm + w0*Yjm0) + (1.0 - upwind) * (w*Yj + w0*Yj0);
                    Hs = upwind * (w*Hjm + w0*Hjm0) + (1.0 - upwind) * (w*Hj + w0*Hj0);
                    
                    g->rho[mc]  = Ys * tmpy;
                    g->rhou[mc] = us * tmpy;  
                    g->rhov[mc] = vs * tmpy + vs * tmpy; 
                    g->rhow[mc] = ws * tmpy; 
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        g->rhoX[nsp][mc] = Xs[nsp] * tmpy; 
                    }
                    g->rhoe[mc] = 0.0 * Hs * tmpy;
                    g->rhoY[mc] = tmpy;
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
            const double oody   = 1.0/elem->dy;
            
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
            
            double oorhoi, ui, vi, wi, Yi, Hi; 
            double oorhoim, uim, vim, wim, Yim, Him;
            double Xi[NSPEC], Xim[NSPEC];
            
            double oorhoj, uj, vj, wj, Yj, Hj; 
            double oorhojm, ujm, vjm, wjm, Yjm, Hjm;
            double Xj[NSPEC], Xjm[NSPEC];
            
            double oorhok, uk, vk, wk, Yk, Hk; 
            double oorhokm, ukm, vkm, wkm, Ykm, Hkm;
            double Xk[NSPEC], Xkm[NSPEC];
            
            double us, vs, ws, Hs, Ys;
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
                        Hi      = oorhoi * (INTERPOL(Sol->rhoe[icp], Sol->rhoe[ic], Sol->rhoe[icm]) + p0_c);
                        
                        oorhoim = 1.0 / INTERPOL(Sol->rhoY[icmm], Sol->rhoY[icm], Sol->rhoY[ic]); 
                        uim     = oorhoim * INTERPOL(Sol->rhou[icmm], Sol->rhou[icm], Sol->rhou[ic]);
                        vim     = oorhoim * INTERPOL(Sol->rhov[icmm], Sol->rhov[icm], Sol->rhov[ic]);
                        wim     = oorhoim * INTERPOL(Sol->rhow[icmm], Sol->rhow[icm], Sol->rhow[ic]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xim[nsp]      = oorhoim * INTERPOL(Sol->rhoX[nsp][icmm], Sol->rhoX[nsp][icm], Sol->rhoX[nsp][ic]);
                        }
                        Yim     = oorhoim * INTERPOL(Sol->rho[icmm], Sol->rho[icm], Sol->rho[ic]);
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
                        Hj = oorhoj * (INTERPOL(Sol->rhoe[jcp], Sol->rhoe[jc], Sol->rhoe[jcm]) + p0_c);
                        
                        oorhojm = 1.0 / INTERPOL(Sol->rhoY[jcmm], Sol->rhoY[jcm], Sol->rhoY[jc]);
                        ujm = oorhojm * INTERPOL(Sol->rhou[jcmm], Sol->rhou[jcm], Sol->rhou[jc]);
                        vjm = oorhojm * INTERPOL(Sol->rhov[jcmm], Sol->rhov[jcm], Sol->rhov[jc]);
                        wjm = oorhojm * INTERPOL(Sol->rhow[jcmm], Sol->rhow[jcm], Sol->rhow[jc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xjm[nsp]      = oorhojm * INTERPOL(Sol->rhoX[nsp][jcmm], Sol->rhoX[nsp][jcm], Sol->rhoX[nsp][jc]);
                        }
                        Yjm = oorhojm * INTERPOL(Sol->rho[jcmm], Sol->rho[jcm], Sol->rho[jc]);
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
                        Hk = oorhok * (INTERPOL(Sol->rhoe[kcp], Sol->rhoe[kc], Sol->rhoe[kcm]) + p0_c);
                        
                        oorhokm = 1.0 / INTERPOL(Sol->rhoY[kcmm], Sol->rhoY[kcm], Sol->rhoY[kc]);
                        ukm = oorhokm * INTERPOL(Sol->rhou[kcmm], Sol->rhou[kcm], Sol->rhou[kc]);
                        vkm = oorhokm * INTERPOL(Sol->rhov[kcmm], Sol->rhov[kcm], Sol->rhov[kc]);
                        wkm = oorhokm * INTERPOL(Sol->rhow[kcmm], Sol->rhow[kcm], Sol->rhow[kc]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Xkm[nsp]      = oorhokm * INTERPOL(Sol->rhoX[nsp][kcmm], Sol->rhoX[nsp][kcm], Sol->rhoX[nsp][kc]);
                        }
                        Ykm = oorhokm * INTERPOL(Sol->rho[kcmm], Sol->rho[kcm], Sol->rho[kc]);
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
					
                    coeff = Gammainv * 0.5 * (Sol->rhoY[ic]*Sol->rhoY[ic]/Sol->rho[ic] + Sol->rhoY[icm]*Sol->rhoY[icm]/Sol->rho[icm]);

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
					
                    coeff = Gammainv * 0.5 * (Sol->rhoY[jc]*Sol->rhoY[jc]/Sol->rho[jc] + Sol->rhoY[jcm]*Sol->rhoY[jcm]/Sol->rho[jcm]);

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

/* ========================================================================== */

void update_SI_MIDPT_buoyancy(ConsVars* Sol, 
                              const ConsVars* flux[3], 
                              const MPV* mpv,
                              const ElemSpaceDiscr* elem,
                              const double dt)
{
    /* 
     implicit midpoint for gravity implies that the buoyancy variable
     receives, besides its own advection update, another update 
     contribution due to the vertical advection of the background
     stratification. That is what we do here based on the flux
     determined at the half time level, and done for a full time
     step.
     */
    const int ndim = elem->ndim;
    const double lambda = dt/elem->dy;
    
    switch (ndim) {
        case 1:
            ERROR("Routine not implemented for 1D\n");
            break;

        case 2: {
            
            const int icx = elem->icx;
            const int igx = elem->igx;
            const int ifx = elem->ifx;
            const int icy = elem->icy;
            const int igy = elem->igy;
            const int ify = elem->ify;
                        
            for (int j=igy; j<icy-igy; j++) {
                int ncj = j*icx;
                int nfj = j;
                double S0p = mpv->HydroState_n->S0[j+1];
                double S0m = mpv->HydroState_n->S0[j];
                double S0c = mpv->HydroState->S0[j];
                for (int i=igx; i<icx-igx; i++) {
                    int ncji = ncj+i;
                    int nfji = nfj+i*ify;
                    Sol->rhoX[BUOY][ncji] += -lambda*(flux[1]->rhoY[nfji+1]*(S0p-S0c) + flux[1]->rhoY[nfji]*(S0c-S0m));
                }
            }
            break;
        }
        case 3:
            ERROR("Routine not implemented for 3D\n");
            break;

        default:
            break;
    }
    
    
}


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log:$
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
