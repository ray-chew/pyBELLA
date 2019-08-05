/*******************************************************************************
 File:   explicit.c
 Author: Rupert, Nicola
 Date:   Mon Mar  2 10:52:55 CET 1998
 *******************************************************************************/
#include <assert.h>
#include <string.h>
#include <math.h>

#include "Common.h"
#include "math_own.h"
#include "boundary.h"
#include "recovery.h"
#include "numerical_flux.h"
#include "error.h"
#include "Eos.h"
#include "userdata.h"
#include "thermodynamic.h"
#include "variable.h"
#include "explicit.h"
#include "enumerator.h"
#include "memory.h"

#if (OUTPUT_SUBSTEPS)
#include "io.h"
#endif

static enum Boolean allocated = WRONG;

static int arraysize = 0;

static States* Lefts;

static States* Rights;

static ConsVars* Fluxes;

static void (*rotate[])(ConsVars* Sol, const enum Direction dir) = {NULL, rotate2D, rotate3D};


/* ============================================================ */

void Explicit_malloc(const int size) {
	extern User_Data ud;
	Fluxes = ConsVars_new(size);
    Lefts = States_new(size);
    Rights = States_new(size);
	allocated = CORRECT; 
	arraysize = size;
}



void Explicit_free() {
	extern User_Data ud;
	ConsVars_free(Fluxes); 
    States_free(Lefts);
    States_free(Rights);
	allocated = WRONG;
	arraysize = 0;
}

/* ================================================================================ */

void copy_fluxes1(double *F, 
                  const double *f, 
                  const int kcache,
                  const int nmax, 
                  const int njump, 
                  const ElemSpaceDiscr* elem)
{
    /* here *f points to the very beginning of the flux array */
    
    for (int i=0; i<nmax; i++) {
        int nc = kcache*njump+i;
        int kc = nc/(elem->icy*elem->icx);
        int jc = (nc-kc*(elem->icy*elem->icx))/elem->icx;
        int ic = nc-kc*(elem->icy*elem->icx)-jc*elem->icx;
        int ii = kc*elem->icy*elem->ifx + jc*elem->ifx + ic + 1;
        F[i]   = f[ii];
    }
}

/* ================================================================================ */

void copy_fluxes(double *F, 
                 const double *f, 
                 const int kcache,
                 const int nmax, 
                 const int njump, 
                 const ElemSpaceDiscr* elem)
{
    const int n0 = kcache*njump;
    const int k0 = n0/(elem->icy*elem->icx);
    const int j0 = (n0-k0*(elem->icy*elem->icx))/elem->icx;
    const int i0 = n0-k0*(elem->icy*elem->icx)-j0*elem->icx;
    
    for (int i=0; i<nmax; i++) {
        int m = (i0+i)/elem->icx;
        F[i] = f[i+m];
    }
}

/* ================================================================================ */

void Explicit_step_and_flux(
							ConsVars* Sol,
							ConsVars* flux,
							const double lambda, 
							const int n, 
							const int SplitStep,
                            const int RK_stage,
                            const enum FluxesFrom adv_fluxes_from, 
                            const enum MUSCL_ON_OFF muscl_on_off) {
	
    /* TODO: Control lambda through  dt  and  ud.tips.update_frac[]
     */
    
	extern User_Data ud;
	extern ElemSpaceDiscr* elem;
	extern ConsVars* dSol;
	extern States* Solk; 
	
	const int ncache = ud.ncache;
	const int njump  = ncache - 2*elem->igx;
    	
	ConsVars pdSol, ppdSol, pFluxes, pflux, ppflux;
    
	double flux_weight_old, flux_weight_new;
    double flux_rhoY_weight_old, flux_rhoY_weight_new;
    	
	int icx, kcache, i, nmax, nsp; 
	int count;
	    
    if (lambda == 0.0) {
        return;
    }
        
    assert(allocated == CORRECT);
	assert(arraysize >= ncache);
    
    icx = elem->icx;
        
    /* bring dummy cells in the current space direction up to date  */
    Bound(Sol, lambda, n, SplitStep);
        
    States_setp(Solk, Sol, 0);
	
	ConsVars_setp(&pdSol, dSol, 0);
	ConsVars_setp(&pflux, flux, 0);
    
	count = 0;
	for (kcache = 0; kcache * njump < n - elem->igx; kcache++) {
		
		nmax = MIN_own(ncache, n - kcache * njump);
		const enum Boolean last = ((kcache + 1) * njump < n - elem->igx) ? WRONG : CORRECT;
		
        if (adv_fluxes_from == FLUX_EXTERNAL) {
            copy_fluxes1(Fluxes->rhoY, flux->rhoY, kcache, nmax, njump, elem);        
        }
        
		/* flux computation*/
        recovery(Lefts, Rights, Solk, Fluxes, lambda, nmax, adv_fluxes_from, muscl_on_off);
        check_flux_bcs(Lefts, Rights, nmax, kcache, njump, elem, SplitStep);
                    
        hllestar(Fluxes, Lefts, Rights, Solk, lambda, nmax, adv_fluxes_from);
		
		/* time updates for conservative variables */
		ConsVars_setp(&ppdSol, &pdSol, elem->igx);
		ConsVars_setp(&ppflux, &pflux, elem->igx);
		ConsVars_setp(&pFluxes, Fluxes, elem->igx-1);
        count  += elem->igx;
		
		/* flux_weight = 0.5; */
		flux_weight_old = ud.tips.flux_frac[RK_stage][0];
		flux_weight_new = ud.tips.flux_frac[RK_stage][1];
        if (adv_fluxes_from == FLUX_EXTERNAL) {
            flux_rhoY_weight_old = 1.0;
            flux_rhoY_weight_new = 0.0;
        } else {
            flux_rhoY_weight_old = ud.tips.flux_frac[RK_stage][0];
            flux_rhoY_weight_new = ud.tips.flux_frac[RK_stage][1];            
        }
        
		for(i = elem->igx; i < nmax - elem->igx; i++) { 
						
            *ppdSol.rho  += lambda * (*pFluxes.rho  - pFluxes.rho[1]);  
			*ppdSol.rhou += lambda * (*pFluxes.rhou - pFluxes.rhou[1]);
			*ppdSol.rhov += lambda * (*pFluxes.rhov - pFluxes.rhov[1]);
			*ppdSol.rhow += lambda * (*pFluxes.rhow - pFluxes.rhow[1]);  
			*ppdSol.rhoe += lambda * (*pFluxes.rhoe - pFluxes.rhoe[1]);
			*ppdSol.rhoY += lambda * (*pFluxes.rhoY - pFluxes.rhoY[1]);
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                *ppdSol.rhoX[nsp] += lambda * (*pFluxes.rhoX[nsp] - pFluxes.rhoX[nsp][1]);
            }
                        
            *ppflux.rho  = flux_weight_old * *ppflux.rho  + flux_weight_new * *pFluxes.rho;  
            *ppflux.rhou = flux_weight_old * *ppflux.rhou + flux_weight_new * *pFluxes.rhou;  
            *ppflux.rhov = flux_weight_old * *ppflux.rhov + flux_weight_new * *pFluxes.rhov;  
            *ppflux.rhow = flux_weight_old * *ppflux.rhow + flux_weight_new * *pFluxes.rhow;  
            *ppflux.rhoe = flux_weight_old * *ppflux.rhoe + flux_weight_new * *pFluxes.rhoe;  
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                *ppflux.rhoX[nsp] = flux_weight_old * *ppflux.rhoX[nsp] + flux_weight_new * *pFluxes.rhoX[nsp];
            }
            *ppflux.rhoY = flux_rhoY_weight_old * *ppflux.rhoY + flux_rhoY_weight_new * *pFluxes.rhoY;  
			
			pFluxes.rho++;
			pFluxes.rhou++;
			pFluxes.rhov++;
			pFluxes.rhow++;
			pFluxes.rhoe++;
			pFluxes.rhoY++;
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                pFluxes.rhoX[nsp]++;
            }
			
			count++;
			ConsVars_addp(&ppdSol, 1);
			ConsVars_addp(&ppflux, 1);
			if(count % icx == 0) {
				ConsVars_addp(&ppflux, 1);
				ConsVars_addp(&pflux, 1);
			}
		}
		if(last == WRONG) {
			ConsVars_addp(&pdSol, njump);
			ConsVars_addp(&pflux, njump); 
			States_addp(Solk, njump); 
			count -= elem->igx;
		}
		else { 
            *ppflux.rho  = flux_weight_old * *ppflux.rho  + flux_weight_new * *pFluxes.rho;  
			*ppflux.rhou = flux_weight_old * *ppflux.rhou + flux_weight_new * *pFluxes.rhou;  
			*ppflux.rhov = flux_weight_old * *ppflux.rhov + flux_weight_new * *pFluxes.rhov;  
			*ppflux.rhow = flux_weight_old * *ppflux.rhow + flux_weight_new * *pFluxes.rhow;  
			*ppflux.rhoe = flux_weight_old * *ppflux.rhoe + flux_weight_new * *pFluxes.rhoe;  
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                *ppflux.rhoX[nsp] = flux_weight_old * *ppflux.rhoX[nsp] + flux_weight_new * *pFluxes.rhoX[nsp];
            }
            *ppflux.rhoY = flux_rhoY_weight_old * *ppflux.rhoY + flux_rhoY_weight_new * *pFluxes.rhoY;  
		}
	}
	
    if (ud.advec_time_integrator == STRANG) {
        Explicit_step_update(Sol, n); 
    }
    
    /* bring dummy cells in the current space direction up to date  */
    Bound(Sol, lambda, n, SplitStep);
}

/* ================================================================================ */

void Absorber(
			  ConsVars* Sol,
              const ElemSpaceDiscr* elem,
			  const double time,
			  const double dt) {
	
	extern User_Data ud;
	extern MPV* mpv;
    extern ConsVars* dSol;
	
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
	
    const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;
    
	const double dy_b = 0.0 * (ud.ymax-ud.ymin);
	const double dx_l = (10.0/240.0) * (ud.xmax-ud.xmin);  
	const double dx_r = (10.0/240.0) * (ud.xmax-ud.xmin);
	
    double dy_t = (9.0/27.0) * (ud.ymax-ud.ymin);  
    double om_y = ud.t_ref/600.0;                               
    double om_x = ud.t_ref/120.0;                    
    if (ud.hill_shape == SCHLUTOW) {
        om_x = 0.0; 
        dy_t = (13.5/27.0) * (ud.ymax-ud.ymin);     
        om_y = ud.t_ref/200.0;             
    }
    int k, j, i, nk, njk, nijk;
	
    for(k = igz; k < icz - igz; k++) {nk = k*icy*icx;
        for(j = igy; j < icy - igy; j++) {njk = nk + j*icx;
            double y = elem->y[j];
			
			const double u_outer    = ud.wind_speed + y * ud.wind_shear; /* velo_background(time); */
			const double v_outer    = 0.0;
			const double w_outer    = 0.0;
			
			
            for(i = igx; i < icx - igx; i++) {nijk = njk + i;
				
                const double x     = elem->x[i];
                const double alpha_l = MAX_own(0.0, ((ud.xmin + dx_l) - x) / dx_l);
                const double alpha_r = MAX_own(0.0, (x - (ud.xmax-dx_r)) / dx_r);
                const double alpha_x = om_x * MAX_own(alpha_l, alpha_r);
                const double alpha_b = MAX_own(0.0, ((ud.ymin + dy_b) - y) / dy_b);
                const double alpha_t = MAX_own(0.0, (y - (ud.ymax-dy_t)) / dy_t);
				const double alpha_y = om_y*MAX_own(alpha_b, alpha_t);
				const double alpha   = MAX_own(alpha_x, alpha_y); 
     	        const double decay   = 1.0 / (1.0 + alpha * dt); 
				
				const double S_outer    = 1.0/mpv->HydroState->Y0[j];
				double Sold, Snew, uold, unew, vold, vnew, wold, wnew; 
				
				Sold = Sol->rho[nijk]  / Sol->rhoY[nijk];
				uold    = Sol->rhou[nijk] / Sol->rho[nijk];
				vold    = Sol->rhov[nijk] / Sol->rho[nijk];
				wold    = Sol->rhow[nijk] / Sol->rho[nijk];
				
				Snew = S_outer + decay * (Sold  - S_outer);
				unew    = u_outer    + decay * (uold     - u_outer);
				vnew    = v_outer    + decay * (vold     - v_outer);
				wnew    = w_outer    + decay * (wold     - w_outer);
				
				// Sol->rhoY[nijk] = rhoY_outer + decay * (Sol->rhoY[nijk] - rhoY_outer);  
				Sol->rho[nijk]  = Sol->rhoY[nijk]*Snew;  
				Sol->rhou[nijk] = Sol->rho[nijk] * unew;  
				Sol->rhov[nijk] = Sol->rho[nijk] * vnew;  
				Sol->rhow[nijk] = Sol->rho[nijk] * wnew;  
				
            }
        }
    }
    Set_Explicit_Boundary_Data(Sol, elem);
}


/* ================================================================================ */

void Explicit_step_update(
						  ConsVars* Sol, 
						  const int n) {
	
    extern User_Data ud;
	extern ElemSpaceDiscr* elem;
	extern ConsVars* dSol;
	ConsVars pSol, pdSol;
	
	int i, nsp;
	
	double rhoe_update_factor = 1.0;
    
    /* update conservative variables */
	
	ConsVars_setp(&pSol, Sol, 0);
	ConsVars_setp(&pdSol, dSol, 0);
	for(i = 0; i < n; i++) {		
		*pSol.rho  +=  *pdSol.rho;
		*pSol.rhou +=  *pdSol.rhou;
		*pSol.rhov +=  *pdSol.rhov;
		*pSol.rhow +=  *pdSol.rhow;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            *pSol.rhoX[nsp] +=  *pdSol.rhoX[nsp];
        }
        
		*pSol.rhoe +=  rhoe_update_factor * *pdSol.rhoe;
		*pSol.rhoY +=  *pdSol.rhoY;

        *pdSol.rho   =  0.0; 
		*pdSol.rhou  =  0.0;
		*pdSol.rhov  =  0.0;
		*pdSol.rhow  =  0.0;
		*pdSol.rhoe  =  0.0;
		*pdSol.rhoY  =  0.0;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            *pdSol.rhoX[nsp] = 0.0;
        }
		
		ConsVars_addp(&pSol, 1);
		ConsVars_addp(&pdSol, 1);         
		
	}
}

/* ================================================================================ */

void advect(
            ConsVars *Sol, 
            ConsVars* flux[3],
            const ConsVars *Sol0, 
            const double dt, 
            const ElemSpaceDiscr* elem,
            const enum FluxesFrom adv_fluxes_from, 
            const enum MUSCL_ON_OFF muscl_on_off, 
            const enum No_of_Strang_Sweeps no_of_sweeps,
            const enum TimeIntegrator advec_time_integrator,
            const int odd)
{
    extern User_Data ud;    
    
    if (advec_time_integrator == EXPL_MIDPT || advec_time_integrator == HEUN) {
        
        printf("\n\n====================================================");
        printf("\nAdvection by explicit midpoint rule, dt = %e", dt);
        printf("\n====================================================\n");
        
        int stage = 0;
        double timestep = ud.tips.update_frac[0] * dt;
        
        for(int Split = 0; Split < elem->ndim; Split++) {
            const double lambda = timestep/elem->dx;
            Explicit_step_and_flux(Sol, flux[Split], lambda, elem->nc, Split, stage, adv_fluxes_from, WITHOUT_MUSCL);                
            (*rotate[elem->ndim - 1])(Sol, FORWARD);
        }
        for(int i_OpSplit = 0; i_OpSplit < elem->ndim; i_OpSplit++) {
            (*rotate[elem->ndim - 1])(Sol, BACKWARD);
        }
        fullD_explicit_updates(Sol, Sol, (const ConsVars**)flux, elem, timestep);
        Set_Explicit_Boundary_Data(Sol, elem);      
        
        stage = 1;
        for(int Split = 0; Split < elem->ndim; Split++) {
            const double lambda = timestep/elem->dx;
            Explicit_step_and_flux(Sol, flux[Split], lambda, elem->nc, Split, stage, adv_fluxes_from, WITHOUT_MUSCL);                
            (*rotate[elem->ndim - 1])(Sol, FORWARD);
        }
        for(int i_OpSplit = 0; i_OpSplit < elem->ndim; i_OpSplit++) {
            (*rotate[elem->ndim - 1])(Sol, BACKWARD);
        }
        fullD_explicit_updates(Sol, Sol, (const ConsVars**)flux, elem, ud.tips.update_frac[1]*dt);
        Set_Explicit_Boundary_Data(Sol, elem);      
    }
    else {  /* Strang splitting == STRANG is the default */ 
        /*
         Simple or Strang splitting for the advection step.
         odd = 0:  even time steps
         odd = 1:  odd  time steps
         used to steer alternating Strang sequences to improve symmetries.
         */
        double time_step = (no_of_sweeps == DOUBLE_STRANG_SWEEP ? 0.5*dt : dt);
        
        printf("\n\n====================================================");
        printf("\nAdvection by Strang splitting, dt = %e", dt);
        printf("\n====================================================\n");
        
        int stage = 0;
        if (odd) {
            for(int Split = 0; Split < elem->ndim; Split++) {
                const double lambda = time_step/elem->dx;
                Explicit_step_and_flux(Sol, flux[Split], lambda, elem->nc, Split, stage, adv_fluxes_from, muscl_on_off);                
                (*rotate[elem->ndim - 1])(Sol, FORWARD);
            }
        } else {
            for(int i_OpSplit = 0; i_OpSplit < elem->ndim; i_OpSplit++) {
                int Split = (elem->ndim - 1) - i_OpSplit;
                (*rotate[elem->ndim - 1])(Sol, BACKWARD);
                const double lambda = time_step/elem->dx;
                Explicit_step_and_flux(Sol, flux[Split], lambda, elem->nc, Split, stage, adv_fluxes_from, muscl_on_off);
            }
        }
        
        if (no_of_sweeps == DOUBLE_STRANG_SWEEP) {
            stage = 1;
            if (odd) {
                for(int i_OpSplit = 0; i_OpSplit < elem->ndim; i_OpSplit++) {
                    int Split = (elem->ndim - 1) - i_OpSplit;
                    (*rotate[elem->ndim - 1])(Sol, BACKWARD);
                    const double lambda = time_step/elem->dx;
                    Explicit_step_and_flux(Sol, flux[Split], lambda, elem->nc, Split, stage, adv_fluxes_from, muscl_on_off);
                }
            } else {
                for(int Split = 0; Split < elem->ndim; Split++) {
                    const double lambda = time_step/elem->dx;
                    Explicit_step_and_flux(Sol, flux[Split], lambda, elem->nc, Split, stage, adv_fluxes_from, muscl_on_off);                
                    (*rotate[elem->ndim - 1])(Sol, FORWARD);
                }            
            }
            
        }
        Set_Explicit_Boundary_Data(Sol, elem);
    }
}

/* ================================================================================ */

void fullD_explicit_updates(ConsVars* Sol, 
                            const ConsVars* Sol0,
                            const ConsVars* flux[3], 
                            const ElemSpaceDiscr* elem, 
                            const double dt) 
{
	
	/*
     complete flux-based recomputation of updates; test, whether
     this creates a problem solving the projection step.
	 */
	extern User_Data ud;

	const int icx = elem->icx;
	const int icy = elem->icy;
	const int icz = elem->icz;
	
	const int ifx = elem->ifx;
	const int ify = elem->ify;
	const int ifz = elem->ifz;
	
	const int igx = elem->igx;
	const int igy = elem->igy;
	const int igz = elem->igz;
    
	double lambda_x = dt / elem->dx;
	double lambda_y = dt / elem->dy;
	double lambda_z = dt / elem->dz;
        
    double delta, deltaSol, delmax, ddelmax, ddelmaxu, ddelmaxv, ddelmaxY;
    double deltax, deltay, deltaz;
        
    int i, j, k, nc; 
    int lcx, lfx, lfy, lfz;
    int mcx, mfx, mfy, mfz;
    int ncx, nfx, nfy, nfz;
    int nfxp, nfyp, nfzp;
        
	switch (elem->ndim) {
		case 1:
			ERROR("\n\n1D-case of recompute_updates() not implemented");
			break;
		case 2:
                        
            delmax  = 0.0;
            ddelmax = 0.0;
            ddelmaxu = 0.0;
            ddelmaxv = 0.0;
            ddelmaxY = 0.0;
            
			for (j = igy; j < icy-igy; j++) {
				mcx = j*icx; 
				mfx = j*ifx; 
				mfy = j;
                                
                for (i = igx; i < icx-igx; i++) {
					ncx  = mcx + i; 
					nfx  = mfx + i;
					nfxp = nfx + 1; 
					nfy  = mfy + i*ify;
					nfyp = nfy + 1;
                    
                    nc = ncx;
                    
                    deltax        = - lambda_x * (flux[0]->rho[nfxp]   - flux[0]->rho[nfx]);
					deltay        = - lambda_y * (flux[1]->rho[nfyp]   - flux[1]->rho[nfy]);
					delta         = deltax + deltay;
                    deltaSol      = Sol->rho[nc] - Sol0->rho[nc];
                    ddelmax       = MAX_own(ddelmax, fabs(delta-deltaSol));
                    Sol->rho[nc]  = Sol0->rho[nc] + delta;
                    
                    deltax        = - lambda_x * (flux[0]->rhou[nfxp]   - flux[0]->rhou[nfx]);
					deltay        = - lambda_y * (flux[1]->rhou[nfyp]   - flux[1]->rhou[nfy]);
					delta         = deltax + deltay;
                    deltaSol      = Sol->rhou[nc] - Sol0->rhou[nc];
                    ddelmaxu       = MAX_own(ddelmaxu, fabs(delta-deltaSol));
                    Sol->rhou[nc] = Sol0->rhou[nc] + delta;

                    deltax        = - lambda_x * (flux[0]->rhov[nfxp]   - flux[0]->rhov[nfx]);
					deltay        = - lambda_y * (flux[1]->rhov[nfyp]   - flux[1]->rhov[nfy]);
					delta         = deltax + deltay;
                    deltaSol      = Sol->rhov[nc] - Sol0->rhov[nc];
                    ddelmaxv       = MAX_own(ddelmaxv, fabs(delta-deltaSol));
                    Sol->rhov[nc] = Sol0->rhov[nc] + delta;
                    
                    deltax        = - lambda_x * (flux[0]->rhow[nfxp]   - flux[0]->rhow[nfx]);
					deltay        = - lambda_y * (flux[1]->rhow[nfyp]   - flux[1]->rhow[nfy]);
					delta         = deltax + deltay;
                    Sol->rhow[nc] = Sol0->rhow[nc] + delta;

                    deltax        = - lambda_x * (flux[0]->rhoe[nfxp]   - flux[0]->rhoe[nfx]);
					deltay        = - lambda_y * (flux[1]->rhoe[nfyp]   - flux[1]->rhoe[nfy]);
					delta         = deltax + deltay;
                    Sol->rhoe[nc] = Sol0->rhoe[nc] + delta;

                    for (int ispec = 0; ispec<ud.nspec; ispec++) {
                        deltax               = - lambda_x * (flux[0]->rhoX[ispec][nfxp]   - flux[0]->rhoX[ispec][nfx]);
                        deltay               = - lambda_y * (flux[1]->rhoX[ispec][nfyp]   - flux[1]->rhoX[ispec][nfy]);
                        delta                = deltax + deltay;
                        Sol->rhoX[ispec][nc] = Sol0->rhoX[ispec][nc] + delta;
                    }
                    
					deltax        = - lambda_x * (flux[0]->rhoY[nfxp]   - flux[0]->rhoY[nfx]);
					deltay        = - lambda_y * (flux[1]->rhoY[nfyp]   - flux[1]->rhoY[nfy]);
					delta         = deltax + deltay;
                    deltaSol      = Sol->rhoY[nc] - Sol0->rhoY[nc];
                    delmax        = MAX_own(delmax, fabs(deltaSol));
                    ddelmaxY      = MAX_own(ddelmaxY, fabs(delta-deltaSol));
                    Sol->rhoY[nc] = Sol0->rhoY[nc] + delta;
                }
			}
            
            /* 
            printf("ddelmax, -u, -v, -Y = %e, %e, %e, %e\n", ddelmax, ddelmaxu, ddelmaxv, ddelmaxY); 
             */

			break;
		case 3:
                        
            delmax  = 0.0;
            ddelmax = 0.0;
            
			for (k = igz; k < icz-igz; k++) {
                lcx = k*icy*icx;
                
                lfx = k*icy*ifx;
                lfy = k*ify;
                lfz = k;
                
                for (j = igy; j < icy-igy; j++) {
                    mcx = lcx + j*icx; 
                    
                    mfx = lfx + j*ifx; 
                    mfy = lfy + j;
                    mfz = lfz + j*icx*ifz;
                    
                    for (i = igx; i < icx-igx; i++) {
                        ncx  = mcx + i; 
                        
                        nfx  = mfx + i;
                        nfy  = mfy + i*icz*ify;
                        nfz  = mfz + i*ifz;

                        nfxp = nfx + 1; 
                        nfyp = nfy + 1;
                        nfzp = nfz + 1;
                        
                        nc = ncx;
                        
                        deltax        = - lambda_x * (flux[0]->rho[nfxp]   - flux[0]->rho[nfx]);
                        deltay        = - lambda_y * (flux[1]->rho[nfyp]   - flux[1]->rho[nfy]);
                        deltaz        = - lambda_z * (flux[2]->rho[nfzp]   - flux[2]->rho[nfz]);
                        delta         = deltax + deltay + deltaz;
                        Sol->rho[nc]  = Sol0->rho[nc] + delta;
                        
                        deltax        = - lambda_x * (flux[0]->rhou[nfxp]   - flux[0]->rhou[nfx]);
                        deltay        = - lambda_y * (flux[1]->rhou[nfyp]   - flux[1]->rhou[nfy]);
                        deltaz        = - lambda_z * (flux[2]->rhou[nfzp]   - flux[2]->rhou[nfz]);
                        delta         = deltax + deltay + deltaz;
                        Sol->rhou[nc] = Sol0->rhou[nc] + delta;
                        
                        deltax        = - lambda_x * (flux[0]->rhov[nfxp]   - flux[0]->rhov[nfx]);
                        deltay        = - lambda_y * (flux[1]->rhov[nfyp]   - flux[1]->rhov[nfy]);
                        deltaz        = - lambda_z * (flux[2]->rhov[nfzp]   - flux[2]->rhov[nfz]);
                        delta         = deltax + deltay + deltaz;
                        Sol->rhov[nc] = Sol0->rhov[nc] + delta;
                        
                        deltax        = - lambda_x * (flux[0]->rhow[nfxp]   - flux[0]->rhow[nfx]);
                        deltay        = - lambda_y * (flux[1]->rhow[nfyp]   - flux[1]->rhow[nfy]);
                        deltaz        = - lambda_z * (flux[2]->rhow[nfzp]   - flux[2]->rhow[nfz]);
                        delta         = deltax + deltay + deltaz;
                        Sol->rhow[nc] = Sol0->rhow[nc] + delta;
                        
                        deltax        = - lambda_x * (flux[0]->rhoe[nfxp]   - flux[0]->rhoe[nfx]);
                        deltay        = - lambda_y * (flux[1]->rhoe[nfyp]   - flux[1]->rhoe[nfy]);
                        deltaz        = - lambda_z * (flux[2]->rhoe[nfzp]   - flux[2]->rhoe[nfz]);
                        delta         = deltax + deltay + deltaz;
                        Sol->rhoe[nc] = Sol0->rhoe[nc] + delta;
                        
                        for (int ispec=0; ispec<ud.nspec; ispec++) {
                            deltax               = - lambda_x * (flux[0]->rhoX[ispec][nfxp]   - flux[0]->rhoX[ispec][nfx]);
                            deltay               = - lambda_y * (flux[1]->rhoX[ispec][nfyp]   - flux[1]->rhoX[ispec][nfy]);
                            deltaz               = - lambda_z * (flux[2]->rhoX[ispec][nfzp]   - flux[2]->rhoX[ispec][nfz]);
                            delta                = deltax + deltay + deltaz;
                            Sol->rhoX[ispec][nc] = Sol0->rhoX[ispec][nc] + delta;
                        }
                                                
                        deltax        = - lambda_x * (flux[0]->rhoY[nfxp]   - flux[0]->rhoY[nfx]);
                        deltay        = - lambda_y * (flux[1]->rhoY[nfyp]   - flux[1]->rhoY[nfy]);
                        deltaz        = - lambda_z * (flux[2]->rhoY[nfzp]   - flux[2]->rhoY[nfz]);
                        delta         = deltax + deltay + deltaz;
                        deltaSol      = Sol->rhoY[nc] - Sol0->rhoY[nc];
                        delmax  = MAX_own(delmax,fabs(deltaSol));
                        ddelmax = MAX_own(ddelmax, fabs(delta-deltaSol));
                        Sol->rhoY[nc] = Sol0->rhoY[nc] + delta;
                    }
                }			
			}	
            
            break;
		default:
			ERROR("\n\nwe are not doing string theory");;
	}

#if OUTPUT_SUBSTEPS /* 4 */
    extern NodeSpaceDiscr* node;
    extern int step;
    if (step >= OUTPUT_SUBSTEPS - 1) {
        putout(Sol, ud.file_name, "Sol", elem, node, 1);
    }
#endif

}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: explicit.c,v $
 Revision 1.1  1998/03/07 09:56:45  nicola
 Added flux computation and multiple pressure variables.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/

