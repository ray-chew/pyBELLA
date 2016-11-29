/*******************************************************************************
 File:   explicit.c
 Author: Rupert, Nicola
 Date:   Mon Mar  2 10:52:55 CET 1998
 *******************************************************************************/
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
#include "io.h"
#include "set_ghostcells_p.h"
#include <assert.h>
#include <string.h>
#include <math.h>

static enum Boolean allocated = WRONG;

static int arraysize = 0;

static States* Lefts;

static States* Rights;

static ConsVars* Fluxes;

static double* gravity_source;

/* ============================================================ */

void Explicit_malloc(const int size) {
	extern User_Data ud;
	Fluxes = ConsVars_new(size);
    Lefts = States_new(size);
    Rights = States_new(size);
    gravity_source = (double*)malloc(size * sizeof(double));
	allocated = CORRECT; 
	arraysize = size;
}


/* ================================================================================ */

void Explicit_free() {
	extern User_Data ud;
	ConsVars_free(Fluxes); 
    States_free(Lefts);
    States_free(Rights);
    free(gravity_source);
	allocated = WRONG;
	arraysize = 0;
}

/* ================================================================================ */

void Explicit_step_and_flux(
							ConsVars* Sol,
							ConsVars* flux,
                            double* buoyS,
                            VectorField* buoy,
                            double* dp2,
							const States* HydroState,
							const double lambda, 
							const int n, 
							const int SplitStep,
                            const int RK_stage,
                            const int implicit) {
	
    /* TODO: Can I get away without ever computing the theta-perturbation evolution,
        just modifying the momentum balance to include the semi-implicit effects?
     */
    
	extern User_Data ud;
	extern ElemSpaceDiscr* elem;
	extern ConsVars* dSol;
	extern States* Solk; 
    extern double *W0, *W1, *W2, *Sbg;
	
	const double gravity_strength = ud.gravity_strength[SplitStep];
    const double dh  = elem->dx;
    const double g   = ud.gravity_strength[SplitStep];
    const double Msq = ud.Msq;
    const double dt  = lambda * dh;
	const int ncache = ud.ncache;
	const int njump  = ncache - 2*elem->igx;
	
	ConsVars pdSol, ppdSol, pFluxes, pflux, ppflux;
    double *pbuoy, *ppbuoy, *pbuoyS, *ppbuoyS;
    double *pdp2, *pS, *pSbg;
    
	double flux_weight_old, flux_weight_new;
    	
	int icx, kcache, i, nmax, nsp; 
	int count;
	    
    if (lambda == 0.0) {
        return;
    }
    
    assert(allocated == CORRECT);
	assert(arraysize >= ncache);
    
    icx = elem->icx;
    
    /*
    Asymmetry of wave representation in the advected case: 
     1) try extrapolating pressure to the quarter and three quarter time levels.
     2) try working with old time level pressure and buoyancy
    */
    
    /* bring dummy cells in the current space direction up to date  */
    Bound(Sol, HydroState, lambda, n, SplitStep);

    double *p2_store = W0;
    double *S_ave = W1;

    assert(elem->ndim == 2); /* lateral averaging for 3D not yet implemented */
    
    for (int ic=0; ic<elem->nc; ic++) {
        p2_store[ic] = Sol->rhoZ[PRES][ic];
    }
    
    for(int j=1; j<elem->icy-1; j++) {
        int njk = j*elem->icx;
        for(int i=0; i<elem->icx; i++) {
            int nijk  = njk + i;
            int nijkp = nijk + elem->icx;
            int nijkm = nijk - elem->icx;
            /* selected weights should correspond to weights in the cell-centered Laplacian */
            Sol->rhoZ[PRES][nijk]  = (ud.latw[0]*p2_store[nijkp] + ud.latw[1]*p2_store[nijk] + ud.latw[2]*p2_store[nijkm]);
            S_ave[nijk] = Sol->rho[nijk]/Sol->rhoY[nijk];
        }
    }
    
    States_setp(Solk, Sol, 0);
    nmax = MIN_own(ncache, n);
    States_HydroState(Solk, HydroState, elem, 0, nmax, 0, SplitStep);
	
	ConsVars_setp(&pdSol, dSol, 0);
	ConsVars_setp(&pflux, flux, 0);
	pbuoy   = (SplitStep == 0 ? buoy->x : (SplitStep == 1 ? buoy->y : buoy->z));
	pdp2    = dp2;
    pbuoyS  = buoyS;
    pS   = S_ave;
    pSbg = Sbg;
    
	count = 0;
	for (kcache = 0; kcache * njump < n - elem->igx; kcache++) {
		
		nmax = MIN_own(ncache, n - kcache * njump);
		const enum Boolean last = ((kcache + 1) * njump < n - elem->igx) ? WRONG : CORRECT;
		                
		/* flux computation*/
        recovery_gravity(Lefts, Rights, gravity_source, pbuoy, pS, pSbg, gravity_strength, Solk, Solk->Y, Solk->rhoZ[PRES], dp2, lambda, nmax, RK_stage, implicit);
        check_flux_bcs(Lefts, Rights, nmax, kcache, njump, elem, SplitStep);
                    
        hllestar(Fluxes, Lefts, Rights, Solk, lambda, nmax);
		
		/* time updates for conservative variables */
		ConsVars_setp(&ppdSol, &pdSol, elem->igx);
		ConsVars_setp(&ppflux, &pflux, elem->igx);
        ppbuoyS = &pbuoyS[elem->igx];
        ppbuoy  = &pbuoy[elem->igx];
		count  += elem->igx;
		ConsVars_setp(&pFluxes, Fluxes, elem->igx-1);
		
		/* flux_weight = 0.5; */
		flux_weight_old = ud.tips.flux_frac[RK_stage][0];
		flux_weight_new = ud.tips.flux_frac[RK_stage][1];
				
		for(i = elem->igx; i < nmax - elem->igx; i++) { 
						
			/* pressure gradient and gravity source terms */
#ifdef GRAVITY_IMPLICIT_1
            double dSbgdy     = (Solk->S0[i+1] - Solk->S0[i-1]) / (2.0*dh);
            double Nsqsc      = - implicit * dt*dt * (g/Msq) * Solk->Y[i] * dSbgdy;
            double ooopNsqsc  = 1.0 / (1.0 + Nsqsc); 
            /* gravity_source[i] = ooopNsqsc * (gravity_source[i] - Nsqsc * Solk->rhou[i]) / lambda; */
            gravity_source[i] = ooopNsqsc * (gravity_source[i] - 2.0 * Nsqsc * Solk->rhou[i]) / lambda; 
#endif
            *ppbuoy             = flux_weight_old * *ppbuoy + flux_weight_new * gravity_source[i];
			*ppdSol.rhou       += lambda * gravity_source[i];

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
            *ppflux.rhoY = flux_weight_old * *ppflux.rhoY + flux_weight_new * *pFluxes.rhoY;  
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                *ppflux.rhoX[nsp] = flux_weight_old * *ppflux.rhoX[nsp] + flux_weight_new * *pFluxes.rhoX[nsp];
            }
			
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
            ppbuoyS++;
            ppbuoy++;
			if(count % icx == 0) {
				ConsVars_addp(&ppflux, 1);
				ConsVars_addp(&pflux, 1);
			}
		}
		if(last == WRONG) {
			ConsVars_addp(&pdSol, njump);
			ConsVars_addp(&pflux, njump); 
            pbuoyS+=njump;
            pbuoy+=njump;
            pS+=njump;
            pSbg+=njump;
            pdp2+=njump;
			States_addp(Solk, njump); 
			States_HydroState(Solk, HydroState, elem, 0, nmax, (kcache+1)*njump, SplitStep);
			count -= elem->igx;
		}
		else { 
            *ppflux.rho  = flux_weight_old * *ppflux.rho  + flux_weight_new * *pFluxes.rho;  
			*ppflux.rhou = flux_weight_old * *ppflux.rhou + flux_weight_new * *pFluxes.rhou;  
			*ppflux.rhov = flux_weight_old * *ppflux.rhov + flux_weight_new * *pFluxes.rhov;  
			*ppflux.rhow = flux_weight_old * *ppflux.rhow + flux_weight_new * *pFluxes.rhow;  
			*ppflux.rhoe = flux_weight_old * *ppflux.rhoe + flux_weight_new * *pFluxes.rhoe;  
			*ppflux.rhoY = flux_weight_old * *ppflux.rhoY + flux_weight_new * *pFluxes.rhoY;  
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                *ppflux.rhoX[nsp] = flux_weight_old * *ppflux.rhoX[nsp] + flux_weight_new * *pFluxes.rhoX[nsp];
            }
		}
	}
	
    if (ud.time_integrator == OP_SPLIT || ud.time_integrator == OP_SPLIT_MD_UPDATE) {
        Explicit_step_update(Sol, n); 
    }    

    for (int ic=0; ic<elem->nc; ic++) {
        Sol->rhoZ[PRES][ic] = p2_store[ic];
    }
}

/* ================================================================================ */

void Absorber(
			  ConsVars* Sol,
			  double time,
			  double dt) {
	
	extern User_Data ud;
    extern ElemSpaceDiscr* elem;
	extern MPV* mpv;
    extern ConsVars* dSol;
	
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
	
    const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;
    
	const double dy_b = 0.0 * (ud.ymax-ud.ymin);
	const double dy_t = (9.0/27.0) * (ud.ymax-ud.ymin);  
	const double om_y = ud.t_ref/600.0;                   
	const double dx_l = (10.0/240.0) * (ud.xmax-ud.xmin);  
	const double dx_r = (10.0/240.0) * (ud.xmax-ud.xmin);
	const double om_x = ud.t_ref/120.0;                    
	
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

void fullD_explicit_updates(ConsVars* Sol, 
                            ConsVars* Sol0,
                            ConsVars* flux[3], 
                            double* buoyS,
                            VectorField* buoy,
                            const ElemSpaceDiscr* elem, 
                            const double dt,
                            const int RK_stage) 
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
    
	double lambda_x = ud.tips.update_frac[RK_stage] * dt / elem->dx;
	double lambda_y = ud.tips.update_frac[RK_stage] * dt / elem->dy;
	double lambda_z = ud.tips.update_frac[RK_stage] * dt / elem->dz;
        
    double delta, deltaSol, delmax, ddelmax, ddelmaxu, ddelmaxv, ddelmaxY;
    double deltax, deltay, deltaz;
        
    int i, j, k, nc; 
    int lcx, lcy, lfx, lfy, lfz;
    int mcx, mcy, mfx, mfy, mfz;
    int ncx, ncy, nfx, nfy, nfz;
    int nfxp, nfyp, nfzp;
    
    double buoymax = 0.0; 
    
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
				mcy = j; 
				mfx = j*ifx; 
				mfy = j;
                                
                for (i = igx; i < icx-igx; i++) {
					ncx  = mcx + i; 
                    ncy  = mcy + i*icy;
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
                    delta        += lambda_x * buoy->x[ncx];
                    deltaSol      = Sol->rhou[nc] - Sol0->rhou[nc];
                    ddelmaxu       = MAX_own(ddelmaxu, fabs(delta-deltaSol));
                    Sol->rhou[nc] = Sol0->rhou[nc] + delta;

                    deltax        = - lambda_x * (flux[0]->rhov[nfxp]   - flux[0]->rhov[nfx]);
					deltay        = - lambda_y * (flux[1]->rhov[nfyp]   - flux[1]->rhov[nfy]);
					delta         = deltax + deltay;
                    delta        += lambda_y * buoy->y[ncy];
                    deltaSol      = Sol->rhov[nc] - Sol0->rhov[nc];
                    ddelmaxv       = MAX_own(ddelmaxv, fabs(delta-deltaSol));
                    Sol->rhov[nc] = Sol0->rhov[nc] + delta;
                    
                    buoymax = MAX_own(buoymax, fabs(buoy->y[ncy]));

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
            printf("buoymax = %e\n", buoymax); 
             
            memset(buoyS, 0.0, elem->nc*sizeof(double));

			break;
		case 3:
            
            delmax  = 0.0;
            ddelmax = 0.0;
            
			for (k = igz; k < icz-igz; k++) {
                lcx = k*icy*icx;
                lcy = k*icy;
                
                lfx = k*icy*ifx;
                lfy = k*ify;
                lfz = k;
                
                for (j = igy; j < icy-igy; j++) {
                    mcx = lcx + j*icx; 
                    mcy = lcy + j;
                    
                    mfx = lfx + j*ifx; 
                    mfy = lfy + j;
                    mfz = lfz + j*icx*ifz;
                    
                    for (i = igx; i < icx-igx; i++) {
                        ncx  = mcx + i; 
                        ncy  = mcy + i*icz*icy;
                        
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
                        delta        += lambda_y * buoy->y[ncy];
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
                        
                        for (int ispec; ispec<ud.nspec; ispec++) {
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
            
            memset(buoyS, 0.0, elem->nc*sizeof(double));

            break;
		default:
			ERROR("\n\nwe are not doing string theory");;
	}
}

#ifdef GRAVITY_IMPLICIT_2
/*------------------------------------------------------------------------------
 explicit step for the fast linear system
 ------------------------------------------------------------------------------*/

void Explicit_Buoyancy(ConsVars* Sol, 
                       VectorField* buoy, 
                       const MPV* mpv, 
                       const ElemSpaceDiscr* elem, 
                       const NodeSpaceDiscr* node, 
                       const double t, 
                       const double dt,
                       const int implicit)
{
    extern User_Data ud;
    extern Thermodynamic th;
    extern double *W0;
    
    double *vold = W0;
    
    double g        = ud.gravity_strength[1];
    double Msq      = ud.Msq;
    double Gamma    = th.Gamma;
    double Gammainv = th.Gammainv;
    
    /* variant relying on the nodal pressure */
    double* p2 = mpv->p2_nodes;
    
    const int icxe = elem->icx;
    const int icye = elem->icy;
    const int icze = elem->icz;
    
    const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;
    
    const int icxn = node->icx;
    const int icyn = node->icy;
    
    const double dx = elem->dx;
    const double dy = elem->dy;
    
    const int kkmax = (elem->ndim > 2 ? 2 : 1);
    const int jjmax = (elem->ndim > 1 ? 2 : 1);
    const int iimax = (elem->ndim > 0 ? 2 : 1);
    
    const int strkn = icyn*icxn;
    const int strjn = icxn;
    const int strin = 1;
    
    double dpdx_max = 0.0;
                
    assert(elem->ndim == 2);

    for (int k=igz; k<icze-igz; k++) {
        int nck = k*icye*icxe;
        int nnk = k*icyn*icxn;
        
        for (int j=igy; j<icye-igy; j++) {
            int nckj   = nck + j*icxe;
            int nckj_y = nck + j;
            int nnkj   = nnk + j*icxn;
            
            /* scaled Brunt-Väisälä frequency - squared */
            double dSdy      = (mpv->HydroState_n->S0[j+1] - mpv->HydroState_n->S0[1]) / dy;
            double dtsqgdSdz = implicit * dt*dt * (g/Msq) * dSdy;  
            
            for (int i=igx; i<icxe-igx; i++) {
                int nckji   = nckj + i;
                int nckji_y = nckj_y + i*icye;
                int nnkji   = nnkj + i;
                
                int nsw = nckji - icxe - 1;
                int nsc = nckji - icxe;
                int nse = nckji - icxe + 1;

                int ncw = nckji - 1;
                int ncc = nckji;
                int nce = nckji + 1;

                int nnw = nckji + icxe - 1;
                int nnc = nckji + icxe;
                int nne = nckji + icxe + 1;

                double dpdx = 0.0;
                double dpdy = 0.0;
                double dpdy_hy = 0.0;

                double Nsqsc = dtsqgdSdz * Sol->rho[nckji]/Sol->rhoY[nckji];
                
                double Sse, Ssc, Ssw, Sce, Scc, Scw, Sne, Snc, Snw, dpie, dpic, dpiw, dpibg;

                /* note: v_old at entry time level could also be computed from
                 the vertical average of P-fluxes */
                
                /* vertical pressure gradient */
                int cnt = 0;
                for (int kk=0; kk<kkmax; kk++) {
                    for (int ii=0; ii<iimax; ii++) {
                        int nnkji_m = nnkji + kk*strkn + ii*strin;
                        int nnkji_p = nnkji + kk*strkn + ii*strin + strjn;
                        
                        dpdy += p2[nnkji_p] - p2[nnkji_m];
                        cnt++;
                    }
                }
                dpdy /= (cnt*dy);
                
                /* horizontal pressure gradient */
                cnt = 0;
                for (int kk=0; kk<kkmax; kk++) {
                    for (int jj=0; jj<jjmax; jj++) {
                        int nnkji_l = nnkji + kk*strkn + jj*strjn;
                        int nnkji_r = nnkji + kk*strkn + jj*strjn + strin;
                        
                        dpdx += p2[nnkji_r] - p2[nnkji_l];
                        cnt++;
                    }
                }
                dpdx /= (cnt*dx);
                dpdx_max = MAX_own(dpdx_max, fabs(dpdx));
                                
                /* hydrostatic pressure gradient 
                Sse    = (Sol->rhoX[BUOY][nse]/Sol->rho[nse] + mpv->HydroState->S0[j-1]); 
                Ssc    = (Sol->rhoX[BUOY][nsc]/Sol->rho[nsc] + mpv->HydroState->S0[j-1]); 
                Ssw    = (Sol->rhoX[BUOY][nsw]/Sol->rho[nsw] + mpv->HydroState->S0[j-1]); 

                Sce    = (Sol->rhoX[BUOY][nce]/Sol->rho[nce] + mpv->HydroState->S0[j]); 
                Scc    = (Sol->rhoX[BUOY][ncc]/Sol->rho[ncc] + mpv->HydroState->S0[j]); 
                Scw    = (Sol->rhoX[BUOY][ncw]/Sol->rho[ncw] + mpv->HydroState->S0[j]); 
                
                Sne    = (Sol->rhoX[BUOY][nne]/Sol->rho[nne] + mpv->HydroState->S0[j+1]); 
                Snc    = (Sol->rhoX[BUOY][nnc]/Sol->rho[nnc] + mpv->HydroState->S0[j+1]); 
                Snw    = (Sol->rhoX[BUOY][nnw]/Sol->rho[nnw] + mpv->HydroState->S0[j+1]); 

                dpie   = -Gamma*(g/Msq)*0.25*(Sse + 2.0*Sce + Sne);
                dpic   = -Gamma*(g/Msq)*0.25*(Ssc + 2.0*Scc + Snc);
                dpiw   = -Gamma*(g/Msq)*0.25*(Ssw + 2.0*Scw + Snw);
                                
                dpdy_hy      = 0.25*(dpie + 2.0*dpic + dpiw);                 
                */
                
                Sse      = Sol->rhoX[BUOY][nse]/Sol->rho[nse]; 
                Ssc      = Sol->rhoX[BUOY][nsc]/Sol->rho[nsc]; 
                Ssw      = Sol->rhoX[BUOY][nsw]/Sol->rho[nsw]; 
                
                Sce      = Sol->rhoX[BUOY][nce]/Sol->rho[nce]; 
                Scc      = Sol->rhoX[BUOY][ncc]/Sol->rho[ncc]; 
                Scw      = Sol->rhoX[BUOY][ncw]/Sol->rho[ncw]; 
                
                Sne      = Sol->rhoX[BUOY][nne]/Sol->rho[nne]; 
                Snc      = Sol->rhoX[BUOY][nnc]/Sol->rho[nnc]; 
                Snw      = Sol->rhoX[BUOY][nnw]/Sol->rho[nnw]; 
                
                dpie     = 0.25*(Sse + 2.0*Sce + Sne);
                dpic     = 0.25*(Ssc + 2.0*Scc + Snc);
                dpiw     = 0.25*(Ssw + 2.0*Scw + Snw);
                dpibg    = 0.25*(mpv->HydroState->S0[j+1] + 2.0*mpv->HydroState->S0[j] + mpv->HydroState->S0[j-1]);
                
                dpdy_hy  = 0.25*(dpie + dpiw);     
                dpdy_hy += 0.5 *dpic;
                dpdy_hy += dpibg;
                dpdy_hy *= -Gamma*(g/Msq);

                
                Sol->rhou[nckji]       -= dt*Gammainv*Sol->rhoY[nckji]*dpdx;
                vold[nckji]             = Sol->rhov[nckji] / Sol->rho[nckji];
                Sol->rhov[nckji]       += (-dt*Gammainv*Sol->rhoY[nckji]*(dpdy-dpdy_hy) - Sol->rhov[nckji] * Nsqsc) / (1.0 + Nsqsc);                    
                Sol->rhoX[BUOY][nckji] -= dt * Sol->rhov[nckji] * dSdy;  

                buoy->x[nckji]    = 0.0;
                buoy->y[nckji_y]  = 0.0;
            }
        }
    }

    /* advection of the background inverse potential temperature
    for (int k=igz; k<icze-igz; k++) {
        int nck = k*icye*icxe;
        
        for (int j=igy; j<icye-igy; j++) {
            int nckj   = nck + j*icxe;
            
            for (int i=igx; i<icxe-igx; i++) {
                int nckji   = nckj + i;
                
                                
                Sol->rhoX[BUOY][nckji] -= Sol->rho[nckji] * dt * vold[nckji] *  (mpv->HydroState_n->S0[j+1] - mpv->HydroState_n->S0[j]) / dy;
            }
        }
    }
    */
}
#endif


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: explicit.c,v $
 Revision 1.1  1998/03/07 09:56:45  nicola
 Added flux computation and multiple pressure variables.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/

