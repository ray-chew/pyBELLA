/*******************************************************************************
 File:   recovery.c
 Author: Rupert, Nicola
 Date:   
 *******************************************************************************/
#include "Common.h"
#include "variable.h"
#include "math_own.h"
#include "error.h"   
#include "Eos.h"
#include "numerical_flux.h"
#include "userdata.h"
#include "thermodynamic.h"
#include <assert.h>
#include "limiter.h"
#include "recovery.h"
#include "enumerator.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>

static enum Boolean allocated = WRONG;

static int arraysize = 0;

static States* Diffs;

static Characters* Slopes;

static Characters* Ampls;

static void slopes(States* Sol, 
                   const int n);

/* ========================================================================== */

double (*limiter[])(const double a, 
					const double b, 
					const double k) = {
    None,
    MinMod,
    VanLeer,
    VanLeerS,
    Superbee,
    MonotonizedCentral,
    SwebyMunz,
    Rupe,
    NoSlope
};

/* ========================================================================== */

void recovery(States* Lefts,
              States* Rights,
              States* Sol,
              ConsVars* Fluxes,
              const double lambda_input, 
              const int nmax,
              const enum FluxesFrom adv_fluxes_from, 
              const enum MUSCL_ON_OFF muscl_on_off) {
    
    extern User_Data ud;
    extern Thermodynamic th;
    extern ElemSpaceDiscr* elem;
    
    const double gamm   = th.gamm;
    const double lambda = (muscl_on_off == 1 ? lambda_input : 0.0);
    
    int OrderTwo  = ((ud.recovery_order == SECOND) ? 1 : 0);
    
    const double internal_flux = (adv_fluxes_from == FLUX_INTERNAL ? 1.0 : 0.0);
    
    double u[nmax], ul[nmax], ur[nmax];
	int i, nsp;  
    
	assert(allocated == CORRECT); 
	
	assert(arraysize >= nmax);
	
	/* obtain primitive variables for the short States-vector */
	primitives(Sol, 0, nmax);
		
    for( i = 1; i < nmax-1; i++ ) { 
        u[i]  = internal_flux * Sol->u[i] + (1.0-internal_flux) * (0.5*(Fluxes->rhoY[i]+Fluxes->rhoY[i-1])/Sol->rhoY[i]);    
        ul[i] = internal_flux * Sol->u[i] + (1.0-internal_flux) * Fluxes->rhoY[i-1]/(0.5*(Sol->rhoY[i]+Sol->rhoY[i-1]));    
        ur[i] = internal_flux * Sol->u[i] + (1.0-internal_flux) * Fluxes->rhoY[i]/(0.5*(Sol->rhoY[i]+Sol->rhoY[i+1]));    
    }
    
	/* differences of primitive quantities */
	for( i = 0; i < nmax - 1;  i++) {      
		double Yrinv = 1.0/Sol->Y[i+1];
		double Ylinv = 1.0/Sol->Y[i];
		
		Diffs->u[i] =  Sol->u[i+1]*Yrinv - Sol->u[i]*Ylinv;
		Diffs->v[i] =  Sol->v[i+1]*Yrinv - Sol->v[i]*Ylinv;
		Diffs->w[i] =  Sol->w[i+1]*Yrinv - Sol->w[i]*Ylinv;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Diffs->X[nsp][i] =  Sol->X[nsp][i+1]*Yrinv - Sol->X[nsp][i]*Ylinv;
        }
		Diffs->Y[i] =  Yrinv - Ylinv;
    }
    
	/* Projection on right eigenvectors */
	slopes(Sol, nmax);
    
	/* right edge states inside cells = left states at cell interfaces */
	
	/* half-timestep */
	for( i = 1; i < nmax-1; i++ ) { 
#ifdef EGDE_VELOCITIES_IN_MUSCL_STEP
        double uu = ur[i];
#else
        double uu = u[i];
#endif
		Ampls->entro[i]   = 0.5 * Slopes->entro[i] * ( 1 - lambda * uu ); /* entro-entry abused for u */
		Ampls->v[i]       = 0.5 * Slopes->v[i]     * ( 1 - lambda * uu );
		Ampls->w[i]       = 0.5 * Slopes->w[i]     * ( 1 - lambda * uu );
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Ampls->X[nsp][i]     = 0.5 * Slopes->X[nsp][i] * ( 1 - lambda * uu );
        }                                
		Ampls->Y[i]     = 0.5 * Slopes->Y[i] * ( 1 - lambda * uu );
	}
	
	for( i = 1; i < nmax-1; i++ ) {
		double S  = 1.0/Sol->Y[i];
		double Yleft = 1.0 / (S + OrderTwo * Ampls->Y[i]);
		
		Lefts->u[i]   = (Sol->u[i]*S   + OrderTwo * Ampls->entro[i]) * Yleft;
		Lefts->v[i]   = (Sol->v[i]*S   + OrderTwo * Ampls->v[i]) * Yleft;
		Lefts->w[i]   = (Sol->w[i]*S   + OrderTwo * Ampls->w[i]) * Yleft;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Lefts->X[nsp][i]   = (Sol->X[nsp][i]*S   + OrderTwo * Ampls->X[nsp][i]) * Yleft;
        }                                
		Lefts->Y[i]   = Yleft;
	}
    
	/* left edge states inside cells = right states at cell interfaces */
	
	/* half-timestep */  
	for( i = 1; i < nmax-1; i++ ) {
#ifdef EGDE_VELOCITIES_IN_MUSCL_STEP
        double uu = ul[i];
#else
        double uu = u[i];
#endif
		Ampls->entro[i] = -0.5 * Slopes->entro[i] * ( 1 + lambda * uu );
		Ampls->v[i]     = -0.5 * Slopes->v[i]     * ( 1 + lambda * uu );
		Ampls->w[i]     = -0.5 * Slopes->w[i]     * ( 1 + lambda * uu );
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Ampls->X[nsp][i]     = -0.5 * Slopes->X[nsp][i] * ( 1 + lambda * uu );
        }                                
		Ampls->Y[i]     = -0.5 * Slopes->Y[i] * ( 1 + lambda * uu );
	}
	
	for( i = 1; i < nmax-1; i++ ) {
		double S   = 1.0/Sol->Y[i];
		double Yright = 1.0 / (S + OrderTwo * Ampls->Y[i]);
		
		Rights->u[i]   = (Sol->u[i]*S   + OrderTwo * Ampls->entro[i]) * Yright;
		Rights->v[i]   = (Sol->v[i]*S   + OrderTwo * Ampls->v[i]) * Yright;
		Rights->w[i]   = (Sol->w[i]*S   + OrderTwo * Ampls->w[i]) * Yright;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Rights->X[nsp][i]   = (Sol->X[nsp][i]*S   + OrderTwo * Ampls->X[nsp][i]) * Yright;
        }                                
		Rights->Y[i]   = Yright;
	}
	    
    for( i = 0; i < nmax-1;  i++) {        
        Lefts->rhoY[i] = Rights->rhoY[i+1] = 0.5*(Sol->rhoY[i] + Sol->rhoY[i+1]) \
        - OrderTwo * 0.5*lambda*(Sol->u[i+1]*Sol->rhoY[i+1]-Sol->u[i]*Sol->rhoY[i]);
        Lefts->p0[i]   = Rights->p0[i+1]   = pow(Lefts->rhoY[i], gamm);
    }
    
    conservatives_from_uvwYZ(Rights, 1, nmax-1); 
    conservatives_from_uvwYZ(Lefts, 1, nmax-1);
}

/* ========================================================================== */

void recovery_malloc(const int size) {
	Diffs = States_new(size);
	Ampls = Characters_new(size);
	Slopes = Characters_new(size);
	allocated = CORRECT; 
	arraysize = size;
}

/* ========================================================================== */

void recovery_free() {
	States_free(Diffs); 
	Characters_free(Ampls);
	Characters_free(Slopes);
	allocated = WRONG;
	arraysize = 0;
}


/*------------------------------------------------------------------------------
 i) projects state differences between neighboring cells onto the right eigen
 vectors of the local flux jacobian matrix (formally these are amplitudes 
 of the waves active between cells or differences in the (pseudo-) characte
 ristic variables.)
 
 ii) associates with each cell a set of slopes of characteristic variables 
 which are obtained from the results of step i) by a nonlinear limiting 
 procedure
 
 iii) We use Sweby / Munz'  k-limiter (see math.h):
 
 LIMIT( as , bs , sa       , sb       , sab        , k )
 LIMIT( |a|, |b|, sign of a, sign of b, sign of a*b, k )
 
 => MAX_own(0,sab)*sa*MAX_own(MIN_own(k*as,bs),MIN_own(as,k*bs)) 
 ------------------------------------------------------------------------------*/

/* ========================================================================== */

static void slopes(States* Sol, 
                   const int n) {
	
    /* User data */
    extern User_Data ud;
	
    const enum LimiterType limiter_type_velocity = ud.limiter_type_velocity;
    const enum LimiterType limiter_type_scalars  = ud.limiter_type_scalars;
    double kp = ud.kp;
    double kz = ud.kz;
    double kY = ud.kY;  
    
    double aul, aur, avl, avr, awl, awr, aYl, aYr; 
    double aXl[NSPEC], aXr[NSPEC];
	
    double du, dv, dw, dY;
    double dX[NSPEC];
	
    int i, nsp;
	
    for( i = 1; i < n-1; i++ ) {
		
        /* amplitudes of left state difference */
        /* (i) - (i-1) */
        du     = Diffs->u[i-1];
        dv     = Diffs->v[i-1];
        dw     = Diffs->w[i-1];
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            dX[nsp] = Diffs->X[nsp][i-1];
        }                                
        dY     = Diffs->Y[i-1];
		
        aul    = du;
        avl    = dv;
        awl    = dw;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            aXl[nsp] = dX[nsp];
        }                                
        aYl    = dY;
		
        /* characteristic amplitudes of right state difference */
        /* (i+1) - (i) */
        du     = Diffs->u[i];
        dv     = Diffs->v[i];
        dw     = Diffs->w[i];
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            dX[nsp] = Diffs->X[nsp][i];
        }                                
        dY     = Diffs->Y[i];
		
        aur    = du;
        avr    = dv;
        awr    = dw;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            aXr[nsp] = dX[nsp];
        }                                
        aYr    = dY;
		
        /* limited slopes */
        Slopes->entro[i] = (*limiter[limiter_type_velocity])(aul, aur, kp);  /* entro abused for u */
        Slopes->v[i] = (*limiter[limiter_type_velocity])(avl, avr, kz);
        Slopes->w[i] = (*limiter[limiter_type_velocity])(awl, awr, kz);
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Slopes->X[nsp][i] = (*limiter[limiter_type_scalars])(aXl[nsp], aXr[nsp], kz);
        }                                
        Slopes->Y[i] = (*limiter[limiter_type_scalars])(aYl, aYr, kY);
    }
}

/*------------------------------------------------------------------------------
 i) projects state differences between neighboring cells onto the right eigen
 vectors of the local flux jacobian matrix (formally these are amplitudes 
 of the waves active between cells or differences in the (pseudo-) characte
 ristic variables.)
 
 ii) associates with each cell a set of slopes of characteristic variables 
 which are obtained from the results of step i) by a nonlinear limiting 
 procedure
 
 iii) We use Sweby / Munz'  k-limiter (see math.h):
  
 LIMIT( as , bs , sa       , sb       , sab        , k )
 LIMIT( |a|, |b|, sign of a, sign of b, sign of a*b, k )
 
 => MAX_own(0,sab)*sa*MAX_own(MIN_own(k*as,bs),MIN_own(as,k*bs)) 
 ------------------------------------------------------------------------------*/


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
