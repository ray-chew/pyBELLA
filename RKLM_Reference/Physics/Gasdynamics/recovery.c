/*******************************************************************************
 File:   recovery.c
 Author: Rupert, Nicola
 Date:   
 *******************************************************************************/
#include "Common.h"
#include "ProjectionType.h"
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
    
    Sl = Hydros[ic].Yinv[0];
    Sc = Hydros[ic].Yinv[2];
    Sr = Hydros[ic].Yinv[4];
    
    double p   = Hydros[ic].p2[2]*Msq;
    double pi  = pow(p, th.Gamma);
    Hydros[ic].p2c[2] = Hydros[ic].p2[2];
    
#ifdef HYDROSTATIC_RHOY_INTERPOLATION
    double P   = Hydros[ic].rhoY[2]*Msq;
    double piP = pow(P, th.gm1);
    
    for(k=0; k<2; k++) {
        double dhk         = (k-2) * 0.5 * dh; 
        double dpi         = - (th.Gamma*strength) * dhk * 0.5 * ((1+0.5*k)*Sc + (1-0.5*k)*Sl);
        double dp          = pow(pi+dpi,th.Gammainv) - p;
        double dP          = pow(piP+dpi,th.gm1inv) - P;
        Hydros[ic].p2[k]   = Hydros[ic].p2[2] + dp/Msq;
        Hydros[ic].rhoY[k] = Hydros[ic].rhoY[2] + dP;
    }
    for(k=3; k<5; k++) {
        double dhk         = (k-2) * 0.5 * dh; 
        double dpi         = - (th.Gamma*strength) * dhk * 0.5 * ((1+0.5*(4-k))*Sc + (1-0.5*(4-k))*Sr);
        double dp          = pow(pi+dpi,th.Gammainv) - p;
        double dP          = pow(piP+dpi,th.gm1inv) - P;
        Hydros[ic].p2[k]   = Hydros[ic].p2[2] + dp/Msq;
        Hydros[ic].rhoY[k] = Hydros[ic].rhoY[2] + dP;
    }
#else
    for(k=0; k<2; k++) {
        double dhk = (k-2) * 0.5 * dh; 
        double dpi = - (th.Gamma*strength) * dhk * 0.5 * ((1+0.5*k)*Sc + (1-0.5*k)*Sl);
        double dp  = pow(pi+dpi,th.Gammainv) - p;
        Hydros[ic].p2[k] = Hydros[ic].p2[2] + dp/Msq;
    }
    for(k=3; k<5; k++) {
        double dhk  = (k-2) * 0.5 * dh; 
        double dpi = - (th.Gamma*strength) * dhk * 0.5 * ((1+0.5*(4-k))*Sc + (1-0.5*(4-k))*Sr);
        double dp  = pow(pi+dpi,th.Gammainv) - p;
        Hydros[ic].p2[k] = Hydros[ic].p2[2] + dp/Msq;
    }
#endif
    
    /* hydrostatic p-difference from cell-centered theta for stronger feedback on cell-centered momentum */
    {
        double dp, dpi;
        
        dpi = - (th.Gamma*strength) * 0.5 * dh * Sc;
        dp  = pow(pi+dpi,th.Gammainv) - p;
        Hydros[ic].p2c[3] = Hydros[ic].p2c[2] + dp/Msq;
        
        dpi = (th.Gamma*strength) * 0.5 * dh * Sc;
        dp  = pow(pi+dpi,th.Gammainv) - p;
        Hydros[ic].p2c[1] = Hydros[ic].p2c[2] + dp/Msq;                
    }
}


/* ========================================================================== */

void recovery_gravity(
					  States* Lefts,
					  States* Rights,
					  double* gravity_source,
                      double* buoy,
                      double* Yinv_ave,
                      double* Yinvbg,
					  const double strength,
                      States* Sol,
					  double* S2,
					  double* p2,
                      const double* dp2,
					  const double lambda, 
					  const int nmax,
                      const int stage) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	extern ElemSpaceDiscr* elem;
	
	const double Msq   = ud.Msq;
	const double gamm  = th.gamm;
    
	const double dh = elem->dx;
    
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
		Diffs->Y[i] =  Yrinv - Ylinv - (Yinvbg[i+1] - Yinvbg[i]);
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
		Ampls->Y[i]     = 0.5 * (Slopes->Y[i] + 0.5*(Yinvbg[i+1] - Yinvbg[i-1])) * ( 1 - lambda *    u   );
	}
	
	for( i = 1; i < nmax-1; i++ ) {
		double Yinv  = 1.0/Sol->Y[i];
		double Yleft = 1.0 / (Yinv + OrderTwo * Ampls->Y[i]);
		
		Lefts->u[i]   = (Sol->u[i]*Yinv   + OrderTwo * Ampls->entro[i]) * Yleft;
		Lefts->v[i]   = (Sol->v[i]*Yinv   + OrderTwo * Ampls->v[i]) * Yleft;
		Lefts->w[i]   = (Sol->w[i]*Yinv   + OrderTwo * Ampls->w[i]) * Yleft;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Lefts->X[nsp][i]   = (Sol->X[nsp][i]*Yinv   + OrderTwo * Ampls->X[nsp][i]) * Yleft;
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
		Ampls->Y[i]     = -0.5 * (Slopes->Y[i] + 0.5*(Yinvbg[i+1] - Yinvbg[i-1]) )  * ( 1 + lambda *    u   );
	}
	
	for( i = 1; i < nmax-1; i++ ) {
		double Yinv   = 1.0/Sol->Y[i];
		double Yright = 1.0 / (Yinv + OrderTwo * Ampls->Y[i]);
		
		Rights->u[i]   = (Sol->u[i]*Yinv   + OrderTwo * Ampls->entro[i]) * Yright;
		Rights->v[i]   = (Sol->v[i]*Yinv   + OrderTwo * Ampls->v[i]) * Yright;
		Rights->w[i]   = (Sol->w[i]*Yinv   + OrderTwo * Ampls->w[i]) * Yright;
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Rights->X[nsp][i]   = (Sol->X[nsp][i]*Yinv   + OrderTwo * Ampls->X[nsp][i]) * Yright;
        }                                
		Rights->Y[i]   = Yright;
	}
	    
    for( i = 0; i < nmax;  i++) {
        
        int iminus = MAX_own(0,i-1);
        int iplus  = MIN_own(nmax-1,i+1);
        
        Hydros[i].rho[2]  = Sol->rho[i];
        
        Hydros[i].rhoY[0] = Sol->rhoY[iminus];
        Hydros[i].rhoY[1] = 0.5*(Sol->rhoY[iminus] + Sol->rhoY[i]);
        Hydros[i].rhoY[2] = Sol->rhoY[i];
        Hydros[i].rhoY[3] = 0.5*(Sol->rhoY[i] + Sol->rhoY[iplus]);
        Hydros[i].rhoY[4] = Sol->rhoY[iplus];
        
        Hydros[i].Yinv[0] = Yinv_ave[iminus];
        Hydros[i].Yinv[2] = Yinv_ave[i];
        Hydros[i].Yinv[4] = Yinv_ave[iplus];
        Hydros[i].Yinv[1] = 0.5*(Hydros[i].Yinv[0]+Hydros[i].Yinv[2]);
        Hydros[i].Yinv[3] = 0.5*(Hydros[i].Yinv[2]+Hydros[i].Yinv[4]);

        Hydros[i].p2[2]   = Sol->rhoZ[i] / Sol->rho[i];  
        
        HydroStates(Hydros, strength, Msq, i, dh, elem);
        
        Lefts->geopot[i]  = 0.5 * (Sol->geopot[iplus] + Sol->geopot[i]);
        Rights->geopot[i] = 0.5 * (Sol->geopot[i] + Sol->geopot[iminus]);
    }
    
    for( i = 0; i < nmax-1;  i++) {        
        Lefts->rhoY[i] = Rights->rhoY[i+1] = 0.5*(Hydros[i].rhoY[3]+Hydros[i+1].rhoY[1]) 
        - 0.5*lambda*(Sol->u[i+1]*Sol->rhoY[i+1]-Sol->u[i]*Sol->rhoY[i]);
        Lefts->p0[i]   = Rights->p0[i+1]   = pow(Lefts->rhoY[i], gamm);
    }

    conservatives_from_uvwYZ(Rights, 1, nmax-1); 
    conservatives_from_uvwYZ(Lefts, 1, nmax-1);
    
#ifdef BUOYANCY_VIA_OPSPLIT
    for( i = 1; i < nmax; i++ ) {
            gravity_source[i] = 0.0;
    }
#else /* BUOYANCY_VIA_OPSPLIT */
	/* pressure gradient and gravity terms */
    double gss  = 1.0;    
    double gps  = 1.0;
    double gths = THETA_EXP_AT_TIME_LEVEL;

#ifdef EDGE_FOCUSED_BUOYANCY
    for( i = 1; i < nmax; i++ ) {
        double drhou, drhou_left, drhou_right, du;
        
        double u    = 0.5*(Sol->u[i]+Sol->u[i-1]);
        
        {            
            double SlopeY     = (1.0/Sol->Y[i] - 1.0/Sol->Y[i-1]); 
            double uSlopeY    = u*SlopeY;
            double rhoYc      = 0.5 * (Rights->rhoY[i]+Lefts->rhoY[i]);
            double dbuoy_adv  = 0.5 * rhoYc*lambda*dh * (strength/Msq) * uSlopeY * gths;
            
            double dp2hydro_l = 0.5 * ((Hydros[i-1].p2[4]-Hydros[i-1].p2[2]) + (Hydros[i].p2[2]-Hydros[i].p2[0]));
            du                = - 0.5 * (Sol->Z[i] - Sol->Z[i-1] - dp2hydro_l) * gps;
            drhou_right       = du;
            drhou_left        = du;
            drhou             = 0.5*(drhou_right+drhou_left);
            Rights->u[i]     += lambda  * du / Rights->rho[i];
            Lefts->u[i-1]    += lambda  * du / Lefts->rho[i-1];
            Rights->rhou[i]  += lambda  * drhou_right;
            Lefts->rhou[i-1] += lambda  * drhou_left;
            
            // Try building advective update of Yinv already into the hydrostatic computations
            
            gravity_source[i]    = gss * drhou + dbuoy_adv; /* switch grav indep of gss ! */
            gravity_source[i-1] += gss * drhou + dbuoy_adv;            
        }
    }
#else
    for( i = 1; i < nmax; i++ ) {
		double drhou, drhou_left, drhou_right, du;
		
		double u    = Sol->u[i];
		
		{            
            double SlopeY     = 0.5 * (1.0/Sol->Y[i+1] - 1.0/Sol->Y[i-1]); 
            double uSlopeY    = u*SlopeY;
			double rhoYc      = 0.5 * (Rights->rhoY[i]+Lefts->rhoY[i]); 
            double dp2hydro_l = 0.5 * ((Hydros[i-1].p2[4]-Hydros[i-1].p2[2]) + (Hydros[i].p2[2]-Hydros[i].p2[0]));
            du                = - 0.5 * (Sol->Z[i] - Sol->Z[i-1] - dp2hydro_l) * gps;
			drhou_right       = du;
			drhou_left        = du;
			drhou             = 0.5*(drhou_right+drhou_left);
			Rights->u[i]     += lambda  * du / Rights->rho[i];
			Lefts->u[i-1]    += lambda  * du / Lefts->rho[i-1];
			Rights->rhou[i]  += lambda  * drhou_right;
			Lefts->rhou[i-1] += lambda  * drhou_left;

            // Try building advective update of Yinv already into the hydrostatic computations
            
            gravity_source[i]    = gss * drhou + rhoYc*lambda*dh * (strength/Msq) *uSlopeY * gths; /* switch grav indep of gss ! */
			gravity_source[i-1] += gss * drhou;            
		}
	}
#endif
#endif /* BUOYANCY_VIA_OPSPLIT */
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
        
        Hydros[i].Yinv[0] = Sol->rho[iminus]/Sol->rhoY[iminus];
        Hydros[i].Yinv[2] = Sol->rho[i]/Sol->rhoY[i];
        Hydros[i].Yinv[4] = Sol->rho[iplus]/Sol->rhoY[iplus];
        Hydros[i].Yinv[1] = 0.5*(Hydros[i].Yinv[0]+Hydros[i].Yinv[2]);
        Hydros[i].Yinv[3] = 0.5*(Hydros[i].Yinv[2]+Hydros[i].Yinv[4]);
        
        Hydros[i].p2[2]    = Sol->rhoZ[i] / Sol->rho[i];  
        
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
			du          = - 0.5 * (Sol->Z[i] - Sol->Z[i-1] - dp2hydro_left);
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
	
#ifndef LIMITER_FROM_MACROS
    const enum LimiterType limiter_type_velocity = ud.limiter_type_velocity;
    const enum LimiterType limiter_type_scalars  = ud.limiter_type_scalars;
    double kp = ud.kp;
    double kz = ud.kz;
    double kY = ud.kY;  
    double kZ = ud.kZ;
#endif
    
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
		
#ifdef LIMITER_FROM_MACROS
        Slopes->entro[i] = LIMITER_U(aul, aur, kp);  
        Slopes->v[i]     = LIMITER_V(avl, avr, kz);
        Slopes->w[i]     = LIMITER_W(awl, awr, kz);
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Slopes->X[nsp][i] = LIMITER_X(aXl[nsp], aXr[nsp], kz);
        }                                
        Slopes->Y[i]     = LIMITER_Y(aYl, aYr, kY);

#else /* LIMITER_FROM_MACROS */ 
        /* limited slopes */
        Slopes->entro[i] = (*limiter[limiter_type_velocity])(aul, aur, kp);  /* entro abused for u */
        Slopes->v[i] = (*limiter[limiter_type_velocity])(avl, avr, kz);
        Slopes->w[i] = (*limiter[limiter_type_velocity])(awl, awr, kz);
        for (nsp = 0; nsp < ud.nspec; nsp++) {
            Slopes->X[nsp][i] = (*limiter[limiter_type_scalars])(aXl[nsp], aXr[nsp], kz);
        }                                
        Slopes->Y[i] = (*limiter[limiter_type_scalars])(aYl, aYr, kY);
#endif  /* LIMITER_FROM_MACROS */
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
