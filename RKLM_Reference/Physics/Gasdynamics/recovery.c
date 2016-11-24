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

static Hydro* Hydros;

static double* rho_buoy;

static Characters* Ampls;

static void slopes_gravity(
						   States* Sol, 
						   Hydro* Hydros, 
						   const int n);

/* ========================================================================== */

double (*limiter[])(
					const double a, 
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

static void HydroStates(Hydro* Hydros,
                        const double strength,
                        const double Msq,
                        const int ic,
                        const double dh,
                        const ElemSpaceDiscr* elem)
{
    extern Thermodynamic th;
    
    double Sc, Sl, Sr;
    int k;
    
    Sl = Hydros[ic].S[0];
    Sc = Hydros[ic].S[2];
    Sr = Hydros[ic].S[4];
        
    for(k=0; k<2; k++) {
        double dhk       = (k-2) * 0.5 * dh; 
        double dpi       = - (th.Gamma*strength) * dhk * 0.5 * ((1+0.5*k)*Sc + (1-0.5*k)*Sl);
        Hydros[ic].p2[k] = Hydros[ic].p2[2] + dpi/Msq;
    }
    for(k=3; k<5; k++) {
        double dhk       = (k-2) * 0.5 * dh; 
        double dpi       = - (th.Gamma*strength) * dhk * 0.5 * ((1+0.5*(4-k))*Sc + (1-0.5*(4-k))*Sr);
        Hydros[ic].p2[k] = Hydros[ic].p2[2] + dpi/Msq;
    }
}


/* ========================================================================== */

void recovery_gravity(
					  States* Lefts,
					  States* Rights,
					  double* gravity_source,
                      double* buoy,
                      double* S_ave,
                      double* Sbg,
					  const double strength,
                      States* Sol,
					  double* S2,
					  double* p2,
                      const double* dp2,
					  const double lambda, 
					  const int nmax,
                      const int stage,
                      const int implicit) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	extern ElemSpaceDiscr* elem;
	
	const double gamm  = th.gamm;
    
#ifndef GRAVITY_IMPLICIT_2
    const double g     = ud.gravity_strength[1];
    const double Msq   = ud.Msq;
	const double dh    = elem->dx;
    const double dt    = lambda * dh;
#endif
    
    int OrderTwo  = ((ud.recovery_order == SECOND) ? 1 : 0);
        
    double u;
	int i, nsp;  
    
	assert(allocated == CORRECT); 
	
	assert(arraysize >= nmax);
	
	/* obtain primitive variables for the short States-vector */
	primitives(Sol, 0, nmax);
		
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
		Diffs->Y[i] =  Yrinv - Ylinv - ADVECT_THETA_PRIME*(Sbg[i+1] - Sbg[i]);
    }
    
	/* Projection on right eigenvectors */
	slopes_gravity(Sol, Hydros, nmax);
    
	/* right edge states inside cells = left states at cell interfaces */
	
	/* half-timestep */
	for( i = 1; i < nmax-1; i++ ) { 
		u   = Sol->u[i];		
		Ampls->entro[i] = 0.5 * Slopes->entro[i] * ( 1 - lambda *    u   ); /* entro-entry abused for u */
		Ampls->v[i]     = 0.5 * Slopes->v[i]     * ( 1 - lambda *    u   );
		Ampls->w[i]     = 0.5 * Slopes->w[i]     * ( 1 - lambda *    u   );
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Ampls->X[nsp][i]     = 0.5 * Slopes->X[nsp][i] * ( 1 - lambda *    u   );
        }                                
		Ampls->Y[i]     = 0.5 * (Slopes->Y[i] + 0.5*ADVECT_THETA_PRIME*(Sbg[i+1] - Sbg[i-1])) * ( 1 - lambda *    u   );
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
		u   = Sol->u[i];		
		Ampls->entro[i] = -0.5 * Slopes->entro[i] * ( 1 + lambda *    u   );
		Ampls->v[i]     = -0.5 * Slopes->v[i]     * ( 1 + lambda *    u   );
		Ampls->w[i]     = -0.5 * Slopes->w[i]     * ( 1 + lambda *    u   );
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Ampls->X[nsp][i]     = -0.5 * Slopes->X[nsp][i] * ( 1 + lambda *    u   );
        }                                
		Ampls->Y[i]     = -0.5 * (Slopes->Y[i] + 0.5*ADVECT_THETA_PRIME*(Sbg[i+1] - Sbg[i-1]) )  * ( 1 + lambda *    u   );
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

#ifndef GRAVITY_IMPLICIT_2
#ifdef GRAVITY_IMPLICIT_1
    
    /* pressure gradient and gravity terms */
    for( i = 0; i < nmax;  i++) {
        
        int iminus = MAX_own(0,i-1);
        int iplus  = MIN_own(nmax-1,i+1);
        
        Hydros[i].rho[2]  = Sol->rho[i];
        
        Hydros[i].S[0] = S_ave[iminus];
        Hydros[i].S[2] = S_ave[i];
        Hydros[i].S[4] = S_ave[iplus];
        Hydros[i].S[1] = 0.5*(Hydros[i].S[0]+Hydros[i].S[2]);
        Hydros[i].S[3] = 0.5*(Hydros[i].S[2]+Hydros[i].S[4]);

        Hydros[i].Sbg[0] = Sbg[iminus];
        Hydros[i].Sbg[2] = Sbg[i];
        Hydros[i].Sbg[4] = Sbg[iplus];
        Hydros[i].Sbg[1] = 0.5*(Hydros[i].Sbg[0]+Hydros[i].Sbg[2]);
        Hydros[i].Sbg[3] = 0.5*(Hydros[i].Sbg[2]+Hydros[i].Sbg[4]);

        Hydros[i].p2[2]   = Sol->rhoZ[PRES][i];  
        
        HydroStates(Hydros, strength, Msq, i, dh, elem);
        
        Lefts->geopot[i]  = 0.5 * (Sol->geopot[iplus] + Sol->geopot[i]);
        Rights->geopot[i] = 0.5 * (Sol->geopot[i] + Sol->geopot[iminus]);
    }
    
    for( i = 1; i < nmax; i++ ) {
                
        double dbuoy_adv  = 0.0;
        
        double Y          = 0.5*(Sol->rhoY[i]/Sol->rho[i] + Sol->rhoY[i-1]/Sol->rho[i-1]);
        double dSbgdy     = (Sbg[i]-Sbg[i-1])/dh;
        double Nsqsc      = - implicit * dt*dt * (g/Msq) * Y * dSbgdy;
        double ooopNsqsc  = 1.0 / (1.0+Nsqsc);
        
        double rhoY       = 0.5 * (Sol->rhoY[i] + Sol->rhoY[i-1]); 
        double dp2hydro_l = 0.5 * ((Hydros[i-1].p2[4]-Hydros[i-1].p2[2]) + (Hydros[i].p2[2]-Hydros[i].p2[0]));
        
        double rhou       = 0.5 * (Rights->rhou[i] + Lefts->rhou[i-1]);
        double drhou      = - 0.5 * lambda * th.Gammainv * rhoY * (Sol->rhoZ[PRES][i] - Sol->rhoZ[PRES][i-1] - dp2hydro_l);
        double drhou_f    = ooopNsqsc * (drhou - Nsqsc * rhou);
        /* double drhou_f    = ooopNsqsc * (drhou - 2.0 * Nsqsc * rhou); */
                
        Rights->u[i]     += drhou_f / Rights->rho[i];
        Lefts->u[i-1]    += drhou_f / Lefts->rho[i-1];
        Rights->rhou[i]  += drhou_f;
        Lefts->rhou[i-1] += drhou_f;

        Rights->S0[i]  = Hydros[i].Sbg[1];
        Lefts->S0[i-1] = Hydros[i-1].Sbg[3];

        /* weighting of this term for implicit part realized in Explicit_Step_and_Flux() */
        gravity_source[i]    = drhou + dbuoy_adv; 
        gravity_source[i-1] += drhou + dbuoy_adv; 
    }
#else  /* GRAVITY_IMPLICIT_1 */
    /* pressure gradient and gravity terms */
    for( i = 0; i < nmax;  i++) {
        
        int iminus = MAX_own(0,i-1);
        int iplus  = MIN_own(nmax-1,i+1);
        
        Hydros[i].rho[2]  = Sol->rho[i];
        
        Hydros[i].S[0] = S_ave[iminus];
        Hydros[i].S[2] = S_ave[i];
        Hydros[i].S[4] = S_ave[iplus];
        Hydros[i].S[1] = 0.5*(Hydros[i].S[0]+Hydros[i].S[2]);
        Hydros[i].S[3] = 0.5*(Hydros[i].S[2]+Hydros[i].S[4]);

        Hydros[i].p2[2]   = Sol->rhoZ[PRES][i];  
        
        HydroStates(Hydros, strength, Msq, i, dh, elem);
        
        Lefts->geopot[i]  = 0.5 * (Sol->geopot[iplus] + Sol->geopot[i]);
        Rights->geopot[i] = 0.5 * (Sol->geopot[i] + Sol->geopot[iminus]);
    }
        
    for( i = 1; i < nmax; i++ ) {
        
        double u          = 0.5*(Sol->u[i]+Sol->u[i-1]);
        double gps        = th.Gammainv * 0.5 * (Sol->rhoY[i] + Sol->rhoY[i-1]); 
        double SlopeY     = (1.0/Sol->Y[i] - 1.0/Sol->Y[i-1]); 
        double uSlopeY    = u*SlopeY;
        double rhoYc      = 0.5 * (Rights->rhoY[i]+Lefts->rhoY[i]);
        double dbuoy_adv  = 0.5 * rhoYc*lambda*dh * (strength/Msq) * uSlopeY * 0.5;
        
        double dp2hydro_l = 0.5 * ((Hydros[i-1].p2[4]-Hydros[i-1].p2[2]) + (Hydros[i].p2[2]-Hydros[i].p2[0]));
        double drhou      = - 0.5 * (Sol->rhoZ[PRES][i] - Sol->rhoZ[PRES][i-1] - dp2hydro_l) * gps;
        
        Rights->u[i]     += lambda  * drhou / Rights->rho[i];
        Lefts->u[i-1]    += lambda  * drhou / Lefts->rho[i-1];
        Rights->rhou[i]  += lambda  * drhou;
        Lefts->rhou[i-1] += lambda  * drhou;
        gravity_source[i]    = drhou + dbuoy_adv; /* switch grav indep of gss ! */
        gravity_source[i-1] += drhou + dbuoy_adv; 
    }
#endif /* GRAVITY_IMPLICIT_1 */
#endif /* GRAVITY_IMPLICIT_2 */
}

/* ========================================================================== */

void pressure_gradient_and_gravity(
                                   double* gravity_source,
                                   const double strength,
                                   States* Sol,
                                   double* S2,
                                   double* p2,
                                   const double lambda, 
                                   const int nmax) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	extern ElemSpaceDiscr* elem;
	
	const double Msq   = ud.Msq;
    
	double dh = elem->dx;
	double rhoY, u;
	int i;  
    
	assert(allocated == CORRECT); 
	assert(arraysize >= nmax);
	
	/* obtain primitive variables for the short States-vector */
	primitives(Sol, 0, nmax);
	
	for( i = 0; i < nmax;  i++) {
		
	    int iminus = MAX_own(0,i-1);
        int iplus  = MIN_own(nmax-1,i+1);
        
        Hydros[i].rho[2]  = Sol->rho[i];
        
        Hydros[i].rhoY[0] = Sol->rhoY[iminus];
        Hydros[i].rhoY[1] = 0.5*(Sol->rhoY[iminus] + Sol->rhoY[i]);
        Hydros[i].rhoY[2] = Sol->rhoY[i];
        Hydros[i].rhoY[3] = 0.5*(Sol->rhoY[i] + Sol->rhoY[iplus]);
        Hydros[i].rhoY[4] = Sol->rhoY[iplus];
        
        Hydros[i].S[0] = Sol->rho[iminus]/Sol->rhoY[iminus];
        Hydros[i].S[2] = Sol->rho[i]/Sol->rhoY[i];
        Hydros[i].S[4] = Sol->rho[iplus]/Sol->rhoY[iplus];
        Hydros[i].S[1] = 0.5*(Hydros[i].S[0]+Hydros[i].S[2]);
        Hydros[i].S[3] = 0.5*(Hydros[i].S[2]+Hydros[i].S[4]);
        
        Hydros[i].p2[2]    = Sol->rhoZ[PRES][i];  
        
        HydroStates(Hydros, strength, Msq, i, dh, elem);
    }
	
	/* pressure gradient and gravity terms */
	for( i = 1; i < nmax-1; i++ ) {
		double dp2hydro_left = 0.5 * ((Hydros[i-1].p2[4]-Hydros[i-1].p2[2]) + (Hydros[i].p2[2]-Hydros[i].p2[0]));
		double drhou, drhou_left, drhou_right, du;
		
		rhoY = Sol->rhoY[i];
		u    = Sol->u[i];
		
		{            
			/* double SlopeY  = Slopes->Y[i]; */
			double SlopeY  = 0.5 * (1.0/Sol->Y[i+1] - 1.0/Sol->Y[i-1]);
			du          = - 0.5 * (Sol->rhoZ[PRES][i] - Sol->rhoZ[PRES][i-1] - dp2hydro_left);
			drhou_right = du;
			drhou_left  = du;
			drhou       = 0.5*(drhou_right+drhou_left);
			gravity_source[i] = drhou + rhoY*0.5*lambda*dh * (strength/Msq) *u*SlopeY;
			gravity_source[i-1] += drhou;
		}
	}
}

/* ========================================================================== */

void recovery_malloc(const int size) {
	Diffs = States_new(size);
	Ampls = Characters_new(size);
	Slopes = Characters_new(size);
	Hydros = Hydro_new(size);
	rho_buoy = (double*)malloc(size*sizeof(double));
	allocated = CORRECT; 
	arraysize = size;
}

/* ========================================================================== */

void recovery_free() {
	States_free(Diffs); 
	Characters_free(Ampls);
	Characters_free(Slopes);
	Hydro_free(Hydros);
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

static void slopes_gravity(
						   States* Sol, 
						   Hydro* Hydros,
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
