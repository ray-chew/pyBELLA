/*******************************************************************************
 File:   explicit.c
 Author: Rupert, Nicola
 Date:   Mon Mar  2 10:52:55 CET 1998
 *******************************************************************************/
#include "Common.h"
#include "ProjectionType.h"
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
                            VectorField* buoy,
                            double* dp2,
							const States* HydroState,
							const double lambda, 
							const int n, 
							const int SplitStep,
                            const int RK_stage) {
	
	extern User_Data ud;
	extern ElemSpaceDiscr* elem;
	extern ConsVars* dSol;
	extern States* Solk; 
    extern double *W0, *W1, *W2, *Yinvbg;
	
	const double gravity_strength = ud.gravity_strength[SplitStep];
	const int ncache = ud.ncache;
	const int njump = ncache - 2*elem->igx;
	
	ConsVars pdSol, ppdSol, pFluxes, pflux, ppflux;
    double *pbuoy, *ppbuoy;
    double *pdp2, *pYinv, *pYinvbg;
    
	double flux_weight_old, flux_weight_new;
    	
	int icx, kcache, i, nmax, nsp; 
	int count;
	    
    assert(allocated == CORRECT);
	assert(arraysize >= ncache);
    
    icx = elem->icx;
        
    /* bring dummy cells in the current space direction up to date  */
    Bound(Sol, HydroState, lambda, n, SplitStep);

    double *p2_store = W0;
    double *Yinv_ave = W1;

    assert(elem->ndim == 2); /* lateral averaging for 3D not yet implemented */
    
    for (int ic=0; ic<elem->nc; ic++) {
        p2_store[ic] = Sol->rhoZ[ic]/Sol->rho[ic];
    }
    
    for(int j=1; j<elem->icy-1; j++) {
        int njk = j*elem->icx;
        for(int i=0; i<elem->icx; i++) {
            int nijk  = njk + i;
            int nijkp = nijk + elem->icx;
            int nijkm = nijk - elem->icx;
            /* selected weights should correspond to weights in the cell-centered Laplacian */
            Sol->rhoZ[nijk]  = Sol->rho[nijk] * (ud.latw[0]*p2_store[nijkp] + ud.latw[1]*p2_store[nijk] + ud.latw[2]*p2_store[nijkm]);
            Yinv_ave[nijk] = Sol->rho[nijk]/Sol->rhoY[nijk];
        }
    }
    
    States_setp(Solk, Sol, 0);
    nmax = MIN_own(ncache, n);
    States_HydroState(Solk, HydroState, elem, 0, nmax, 0, SplitStep);
	
	ConsVars_setp(&pdSol, dSol, 0);
	ConsVars_setp(&pflux, flux, 0);
	pbuoy = (SplitStep == 0 ? buoy->x : (SplitStep == 1 ? buoy->y : buoy->z));
	pdp2  = dp2;
    pYinv = Yinv_ave;
    pYinvbg  = Yinvbg;
    
	count = 0;
	for (kcache = 0; kcache * njump < n - elem->igx; kcache++) {
		
		nmax = MIN_own(ncache, n - kcache * njump);
		const enum Boolean last = ((kcache + 1) * njump < n - elem->igx) ? WRONG : CORRECT;
		                
		/* flux computation*/
        recovery_gravity(Lefts, Rights, gravity_source, pbuoy, pYinv, pYinvbg, gravity_strength, Solk, Solk->Y, Solk->Z, dp2, lambda, nmax, RK_stage);
        check_flux_bcs(Lefts, Rights, nmax, kcache, njump, elem, SplitStep);
                    
        hllestar(Fluxes, Lefts, Rights, Solk, lambda, nmax);
		
		/* time updates for conservative variables */
		ConsVars_setp(&ppdSol, &pdSol, elem->igx);
		ConsVars_setp(&ppflux, &pflux, elem->igx);
        ppbuoy = &pbuoy[elem->igx];
		count += elem->igx;
		ConsVars_setp(&pFluxes, Fluxes, elem->igx-1);
		
		/* flux_weight = 0.5; */
		flux_weight_old = ud.tips.flux_frac[RK_stage][0];
		flux_weight_new = ud.tips.flux_frac[RK_stage][1];
				
		for(i = elem->igx; i < nmax - elem->igx; i++) { 
			
			*ppflux.rho  = flux_weight_old * *ppflux.rho  + flux_weight_new * *pFluxes.rho;  
			*ppflux.rhou = flux_weight_old * *ppflux.rhou + flux_weight_new * *pFluxes.rhou;  
			*ppflux.rhov = flux_weight_old * *ppflux.rhov + flux_weight_new * *pFluxes.rhov;  
			*ppflux.rhow = flux_weight_old * *ppflux.rhow + flux_weight_new * *pFluxes.rhow;  
			*ppflux.rhoe = flux_weight_old * *ppflux.rhoe + flux_weight_new * *pFluxes.rhoe;  
			*ppflux.rhoY = flux_weight_old * *ppflux.rhoY + flux_weight_new * *pFluxes.rhoY;  
			*ppflux.rhoZ = flux_weight_old * *ppflux.rhoZ + flux_weight_new * *pFluxes.rhoZ;  
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                *ppflux.rhoX[nsp] = flux_weight_old * *ppflux.rhoX[nsp] + flux_weight_new * *pFluxes.rhoX[nsp];
            }
			
			/* pressure gradient and gravity source terms */
            *ppbuoy      = flux_weight_old * *ppbuoy + flux_weight_new * gravity_source[i];
			*ppdSol.rhou += lambda * gravity_source[i];

            *ppdSol.rho  += lambda * (*pFluxes.rho  - pFluxes.rho[1]);  
			*ppdSol.rhou += lambda * (*pFluxes.rhou - pFluxes.rhou[1]);
			*ppdSol.rhov += lambda * (*pFluxes.rhov - pFluxes.rhov[1]);
			*ppdSol.rhow += lambda * (*pFluxes.rhow - pFluxes.rhow[1]);  
			*ppdSol.rhoe += lambda * (*pFluxes.rhoe - pFluxes.rhoe[1]);
			*ppdSol.rhoY += lambda * (*pFluxes.rhoY - pFluxes.rhoY[1]);
			*ppdSol.rhoZ += lambda * (*pFluxes.rhoZ - pFluxes.rhoZ[1]);
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                *ppdSol.rhoX[nsp] += lambda * (*pFluxes.rhoX[nsp] - pFluxes.rhoX[nsp][1]);
            }
			
			pFluxes.rho++;
			pFluxes.rhou++;
			pFluxes.rhov++;
			pFluxes.rhow++;
			pFluxes.rhoe++;
			pFluxes.rhoY++;
			pFluxes.rhoZ++;
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                pFluxes.rhoX[nsp]++;
            }
			
			
			count++;
			ConsVars_addp(&ppdSol, 1);
			ConsVars_addp(&ppflux, 1);
            ppbuoy++;
			if(count % icx == 0) {
				ConsVars_addp(&ppflux, 1);
				ConsVars_addp(&pflux, 1);
			}
		}
		if(last == WRONG) {
			ConsVars_addp(&pdSol, njump);
			ConsVars_addp(&pflux, njump); 
            pbuoy+=njump;
            pYinv+=njump;
            pYinvbg+=njump;
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
			*ppflux.rhoZ = flux_weight_old * *ppflux.rhoZ + flux_weight_new * *pFluxes.rhoZ;  
            for (nsp = 0; nsp < ud.nspec; nsp++) {
                *ppflux.rhoX[nsp] = flux_weight_old * *ppflux.rhoX[nsp] + flux_weight_new * *pFluxes.rhoX[nsp];
            }
		}
	}
	
    if (ud.time_integrator == OP_SPLIT || ud.time_integrator == OP_SPLIT_MD_UPDATE) {
        Explicit_step_update(Sol, n); 
    }    

    for (int ic=0; ic<elem->nc; ic++) {
        Sol->rhoZ[ic] = p2_store[ic]*Sol->rho[ic];
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
				
				const double Yinv_outer    = 1.0/mpv->HydroState->Y0[j];
				double Yinvold, Yinvnew, uold, unew, vold, vnew, wold, wnew; 
				
				Yinvold = Sol->rho[nijk]  / Sol->rhoY[nijk];
				uold    = Sol->rhou[nijk] / Sol->rho[nijk];
				vold    = Sol->rhov[nijk] / Sol->rho[nijk];
				wold    = Sol->rhow[nijk] / Sol->rho[nijk];
				
				Yinvnew = Yinv_outer + decay * (Yinvold  - Yinv_outer);
				unew    = u_outer    + decay * (uold     - u_outer);
				vnew    = v_outer    + decay * (vold     - v_outer);
				wnew    = w_outer    + decay * (wold     - w_outer);
				
				// Sol->rhoY[nijk] = rhoY_outer + decay * (Sol->rhoY[nijk] - rhoY_outer);  
				Sol->rho[nijk]  = Sol->rhoY[nijk]*Yinvnew;  
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
		double rho_old = *pSol.rho;
		
		*pSol.rho  +=  *pdSol.rho;
		*pSol.rhou +=  *pdSol.rhou;
		*pSol.rhov +=  *pdSol.rhov;
		*pSol.rhow +=  *pdSol.rhow;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            *pSol.rhoX[nsp] +=  *pdSol.rhoX[nsp];
        }
        
		*pSol.rhoe +=  rhoe_update_factor * *pdSol.rhoe;
		*pSol.rhoY +=  *pdSol.rhoY;
		*pSol.rhoZ =  (*pSol.rhoZ/rho_old) * *pSol.rho;

        *pdSol.rho   =  0.0; 
		*pdSol.rhou  =  0.0;
		*pdSol.rhov  =  0.0;
		*pdSol.rhow  =  0.0;
		*pdSol.rhoe  =  0.0;
		*pdSol.rhoY  =  0.0;
		*pdSol.rhoZ  =  0.0;
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
                            VectorField* buoy,
                            const ElemSpaceDiscr* elem, 
                            const double dt,
                            const int RK_stage) {
	
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
                    
                    
                    Sol->rhoZ[nc] = (Sol0->rhoZ[nc]/Sol0->rho[nc]) * Sol->rho[nc];

					deltax        = - lambda_x * (flux[0]->rhoY[nfxp]   - flux[0]->rhoY[nfx]);
					deltay        = - lambda_y * (flux[1]->rhoY[nfyp]   - flux[1]->rhoY[nfy]);
					delta         = deltax + deltay;
                    deltaSol      = Sol->rhoY[nc] - Sol0->rhoY[nc];
                    delmax        = MAX_own(delmax,fabs(deltaSol));
                    ddelmaxY      = MAX_own(ddelmaxY, fabs(delta-deltaSol));
                    Sol->rhoY[nc] = Sol0->rhoY[nc] + delta;
                    
                }
			}
            
            /* */
            printf("ddelmax, -u, -v, -Y = %e, %e, %e, %e\n", ddelmax, ddelmaxu, ddelmaxv, ddelmaxY); 
             
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
                        /* delta        += lambda_x * buoy->x[ncx]; */
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
                        /* delta        += lambda_z * buoy->z[ncz];  */
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
                        
                        
                        Sol->rhoZ[nc] = (Sol0->rhoZ[nc]/Sol0->rho[nc]) * Sol->rho[nc];
                        
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
}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: explicit.c,v $
 Revision 1.1  1998/03/07 09:56:45  nicola
 Added flux computation and multiple pressure variables.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/

