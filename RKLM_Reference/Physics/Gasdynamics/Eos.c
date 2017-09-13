/*******************************************************************************
 File:   Eos.c (Equation of state - routines)
 Author: Rupert 
 Date:   ?
 *******************************************************************************/
#include "Common.h"
#include <math.h>
#include "userdata.h"
#include "thermodynamic.h"
#include "variable.h"
#include "mpv.h"
#include "math_own.h"
#include "Eos.h"
#include "set_ghostcells_p.h"


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void primitives(
				States* U, 
				const int nstart, 
				const int nende)
{
	extern User_Data ud;
	extern Thermodynamic th;
	
	const double Msq = ud.compressibility * ud.Msq;
    
	const double gamma = th.gamm;
	double q;
	int i, nsp;
	
	for(i = nstart; i < nende; i++) {
		
		U->u[i] = U->rhou[i] / U->rho[i];
		U->v[i] = U->rhov[i] / U->rho[i];
		U->w[i] = U->rhow[i] / U->rho[i];
		U->Y[i] = U->rhoY[i] / U->rho[i];
        for (nsp=0; nsp<ud.nspec; nsp++) {
            U->X[nsp][i] = U->rhoX[nsp][i] / U->rho[i];
        }
		U->p[i] = pow(U->rhoY[i],gamma);
		q = 0.5 * Msq * (U->u[i] * U->u[i] + U->v[i] * U->v[i] + U->w[i] * U->w[i]);
		U->q[i] = q;
		U->entro[i] = 1.0/U->Y[i];
	}
}


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void conservatives(
				   States* U, 
				   const int nstart, 
				   const int nende)
{
	extern User_Data ud;
	extern Thermodynamic th;
	
	const double Msq = ud.compressibility * ud.Msq;
	const double gm1_inv = 1.0 / th.gm1;
	double q;
	int i, nsp;
	
	for(i = nstart; i < nende; i++) {
		U->rhou[i] = U->u[i] * U->rho[i];
		U->rhov[i] = U->v[i] * U->rho[i];
		U->rhow[i] = U->w[i] * U->rho[i];
		U->rhoY[i] = U->Y[i] * U->rho[i];
        for (nsp=0; nsp<ud.nspec; nsp++) {
            U->rhoX[nsp][i] = U->X[nsp][i] * U->rho[i];
        }
		q = 0.5 * Msq * (U->u[i] * U->u[i] + U->v[i] * U->v[i] + U->w[i] * U->w[i]);
		U->rhoe[i] = gm1_inv * U->p[i] + U->rho[i] * (q + U->geopot[i]);
	}
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void conservatives_from_entro_vars(
								   States* U, 
								   const int nstart, 
								   const int nende)
{
	extern User_Data ud;
	extern Thermodynamic th;
	
	const double Msq = ud.compressibility * ud.Msq;
	const double gm1_inv = 1.0 / th.gm1;
	const double gamm_inv = th.gamminv;
	double q;
	int i, nsp;
	
	for(i = nstart; i < nende; i++) {
		double rho = U->entro[i] * pow(U->p[i], gamm_inv);
		U->rho[i]  = rho;
		U->rhou[i] = U->u[i] * rho;
		U->rhov[i] = U->v[i] * rho;
		U->rhow[i] = U->w[i] * rho;
		U->rhoY[i] = U->Y[i] * rho;
        for (nsp=0; nsp<ud.nspec; nsp++) {
            U->rhoX[nsp][i] = U->X[nsp][i] * U->rho[i];
        }
		q = 0.5 * Msq * (U->u[i] * U->u[i] + U->v[i] * U->v[i] + U->w[i] * U->w[i]);
		U->rhoe[i] = gm1_inv * U->p[i] + rho * (q + U->geopot[i]);
	}
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void conservatives_from_primitives(
								   States* U, 
								   const int nstart, 
								   const int nende)
{
	extern User_Data ud;
	extern Thermodynamic th;
	
	const double Msq = ud.compressibility * ud.Msq;
	const double gm1_inv = 1.0 / th.gm1;
	double q;
	int i, nsp;
	
	for(i = nstart; i < nende; i++) {
		double rho = U->rho[i];
		U->rhou[i] = U->u[i] * rho;
		U->rhov[i] = U->v[i] * rho;
		U->rhow[i] = U->w[i] * rho;
		U->rhoY[i] = U->Y[i] * rho;
        for (nsp=0; nsp<ud.nspec; nsp++) {
            U->rhoX[nsp][i] = U->X[nsp][i] * U->rho[i];
        }
		q = 0.5 * Msq * (U->u[i] * U->u[i] + U->v[i] * U->v[i] + U->w[i] * U->w[i]);
		U->rhoe[i] = gm1_inv * U->p[i] + rho * (q + U->geopot[i]);
	}
}


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void conservatives_from_uvwYZ(
							  States* U, 
							  const int nstart, 
							  const int nende)
{
	extern User_Data ud;
	extern Thermodynamic th;
	
	int i, nsp;
	
	for(i = nstart; i < nende; i++) {
		double rho = U->rhoY[i] / U->Y[i];
		double geopot = U->geopot[i];
        double u = U->u[i];
        double v = U->v[i];
        double w = U->w[i];
		double p;
		
		U->rho[i]  = rho;
		U->rhou[i] = U->u[i] * rho;
		U->rhov[i] = U->v[i] * rho;
		U->rhow[i] = U->w[i] * rho;
		U->rhoY[i] = U->Y[i] * rho;
        for (nsp=0; nsp<ud.nspec; nsp++) {
            U->rhoX[nsp][i] = U->X[nsp][i] * rho;
        }
        
		p = pow(U->rhoY[i] , th.gamminv);  
		U->rhoe[i] = rhoe(rho, u, v, w, p, geopot);
	}
}

/*------------------------------------------------------------------------------
 thermodynamic pressure
 ------------------------------------------------------------------------------*/
void pressure(
			  double *p, 
			  const ConsVars* U, 
			  const int nstart_p, 
			  const int nstart, 
			  const int nende) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	
	int i, ip;
	
	for( ip = nstart_p, i = nstart;  i < nende ; ip++, i++)
    {
		p[ip] = pow(U->rhoY[i],th.gamm);
    }	
}

/*------------------------------------------------------------------------------
 deviation of thermodynamic pressure from background state (dimensional)
 ------------------------------------------------------------------------------*/
void dpress_dim(
                double *dpdim, 
                const ConsVars* U, 
                const MPV *mpv,
                const ElemSpaceDiscr *elem)
{
    extern User_Data ud;
    extern Thermodynamic th;
    
    double p; 
    
    int i, j, k, l, m, n;
    int icx = elem->icx;
    int icy = elem->icy;
    int icz = elem->icz;
    
    for (k=0; k<icz; k++) {
        l = k*icx*icy;
        for (j=0; j<icy; j++) {
            m = l + j*icx;
            for (i=0; i<icx; i++) {
                n = m + i;
                p = pow(th.Gamma * ud.Msq * mpv->p2_cells[n],th.Gammainv);
                dpdim[n] = (p - mpv->HydroState->p0[j]) * ud.p_ref;
            }
        }
    }
}


/*------------------------------------------------------------------------------
 pressure increment computed in the first projection
 ------------------------------------------------------------------------------*/
void dp2_first_projection(
                          double *dp2, 
                          const ConsVars* U, 
                          const MPV *mpv,
                          const ElemSpaceDiscr *elem)
{
    int icx = elem->icx;
    int icy = elem->icy;
    int icz = elem->icz;
    
    for (int k=0; k<icz; k++) {int l = k*icx*icy;
        for (int j=0; j<icy; j++) {int m = l + j*icx;
            for (int i=0; i<icx; i++) {int n = m + i;
                dp2[n] = mpv->dp2_cells[n];
            }
        }
    }
}


/*------------------------------------------------------------------------------
 deviation of Exner pressure from background state values
 ------------------------------------------------------------------------------*/
void dp_exner(
              double *dpi, 
              const ConsVars* U, 
              const MPV *mpv,
              const ElemSpaceDiscr *elem) 
{
    extern User_Data ud;
    extern Thermodynamic th;
    
    int i, j, k, l, m, n;
    int icx = elem->icx;
    int icy = elem->icy;
    int igy = elem->igy;
    int icz = elem->icz;
    
    double p2_o = th.Gammainv * pow(mpv->HydroState->p0[igy],th.Gamma) / ud.Msq;
    
    for (k=0; k<icz; k++) {
        l = k*icx*icy;
        for (j=0; j<icy; j++) {
            m = l + j*icx;
            for (i=0; i<icx; i++) {
                n = m + i;
                dpi[n] = mpv->p2_cells[n] - mpv->HydroState->p20[j] - p2_o;
            }
        }
    }
}


/*------------------------------------------------------------------------------
 inverse of potential temperature
 ------------------------------------------------------------------------------*/
void pot_temp_inv(
				  double *S, 
				  const ConsVars* U, 
				  const int nstart_p, 
				  const int nstart, 
				  const int nende) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	
	int i;
	
	for( i = nstart; i < nende; i++)
    {
        S[i] = U->rho[i] / U->rhoY[i];
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void velox(double *u, const ConsVars* U, const int nstart, const int nende)
{   
	double *prho, *prhou, *pu;
	
	int i;
	
	for( i     = nstart          ,
		prho  = &U->rho[nstart]  ,
		prhou = &U->rhou[nstart] ,
		pu    = &u[nstart]          ;
		i < nende  ;  
		i++     ,
		prho++  ,
		prhou++ , 
		pu++      )
    {
		*pu =  *prhou / *prho;
    }
}


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void veloy(double *v, const ConsVars* U, const int nstart, const int nende)
{   
	double *prho, *prhov, *pv;
	
	int i;
	
	for( i = nstart , prho = &U->rho[nstart], prhov = &U->rhov[nstart], pv = &v[nstart];
		i < nende;  
		i++ ,        prho++,                prhov++,                 pv++)
    {
		*pv =  *prhov / *prho;
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void veloz(double *w, const ConsVars* U, const int nstart, const int nende)
{   
	double *prho, *prhow, *pw;
	
	int i;
	
	for( i = nstart , prho = &U->rho[nstart], prhow = &U->rhow[nstart], pw = &w[nstart];
		i < nende;  
		i++ ,        prho++,                prhow++,                 pw++)
    {
		*pw =  *prhow / *prho;
    }
}


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void fluctuation(double *var, const double *var0, const ElemSpaceDiscr *elem){
    double mean = 0.0;
    int count = 0;
    for (int k=elem->igz; k<elem->icz-elem->igz; k++) {int nk = k*elem->icy*elem->icx;
        for (int j=elem->igy; j<elem->icy-elem->igy; j++) {int njk = nk + j*elem->icy;
            for (int i=elem->igx; i<elem->icx-elem->igx; i++) {int nijk = njk + i;
                mean += var0[nijk];
                count++;
            }
        }
    }
    mean/=count;
    for (int ic=0; ic<elem->nc; ic++) {
        var[ic] = var0[ic] - mean; 
    }
}



/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void delta_rhoY(double *drhoY, const ConsVars* U, const int nstart, const int nende)
{   
    extern MPV* mpv;
    extern ElemSpaceDiscr* elem;
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    const int icz = elem->icz;
        
    for(int k=0; k<icz; k++) {
        int nk = k*icy*icx;
        for (int j=0; j<icy; j++) {
            int nkj = nk + j*icx;
            for (int i=0; i<icx; i++) {
                int nkji = nkj+i;
                drhoY[nkji] =  U->rhoY[nkji] - mpv->HydroState->rho0[j]*mpv->HydroState->Y0[j];
            }
        }
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void MassFrac_Y(double *Y, const ConsVars* U, const int nstart, const int nende)
{   
	double *prho, *prhoY, *pY;  
	int i;
	
	for(i = nstart, 
		prho = &U->rho[nstart],
		prhoY = &U->rhoY[nstart], 
		pY = &Y[nstart];
		i < nende;  
		i++,        
		prho++,                
		prhoY++,                 
		pY++) {
		
		*pY =  *prhoY / *prho;
	}
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void MassFrac_Y_perturbation(double *dY, const ConsVars* U, const int nstart, const int nende)
{   
    extern ElemSpaceDiscr *elem;
    
    int i;
    
    double *prho, *prhoY, *pdY;
	
	for(i = nstart, 
		prho = &U->rho[nstart],
		prhoY = &U->rhoY[nstart], 
		pdY = &dY[nstart];
		i < nende;  
		i++,        
		prho++,                
		prhoY++,                 
		pdY++) {
        
		int k = i/(elem->icy*elem->icx);
        int j = (i-k*elem->icy*elem->icx)/elem->icx;
        double y = elem->y[j];
        double Yb = stratification(y);
		*pdY =  *prhoY / *prho - Yb;
	}
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void MassFrac_q(
                double *q, 
                const ConsVars* U, 
                const int nstart, 
                const int nende,
                const int which_q) {  
    
	int i;
	
	for( i = nstart; i < nende; i++)
    {
		q[i] =  U->rhoX[which_q][i] / U->rho[i];
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void specific_enthalpy(
					   double* h, 
					   const ConsVars* U, 
					   const int nstart, 
					   const int nende) {
	
	double *p;
	int i;
	
	p = h;             /* abuse results-array as intermediate storage */
	pressure(p,U,nstart,nstart,nende);
	
	for(i = nstart; i < nende; i++) {
		h[i] = (U->rhoe[i] + p[i]) / U->rho[i];
	}
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double rhoe(
			const double rho, 
			const double u, 
			const double v, 
			const double w,
			const double p,
			const double geopot) {
	
	extern User_Data ud;
	extern Thermodynamic th;
	
	const double Msq = ud.compressibility * ud.Msq;
	const double gm1inv = th.gm1inv;
	
	return p * gm1inv + 0.5 * Msq * rho * (u * u + v * v + w * w) + rho * geopot;
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void adjust_pi(
			   ConsVars* Sol,
               MPV* mpv,
			   const ConsVars* Sol0,
			   const ElemSpaceDiscr* elem, 
               const double weight) {
	
    /*
     NOTICE: in the fully compressible case I obtain the update for
     Z from the update of rhoY == P; 
     In the psinc case and in case of one of the intermediate models, 
     dp2_elliptic is accounted for as well. In these cases the inter-
     pretation of this quantity is extremely important. Is it a full
     pressure update over one time step or is it half of that or, 
     for the intermediate models, something in between?
     */
    
    
	extern User_Data ud;
	extern Thermodynamic th;
	    
    double alpha = ud.compressibility; 
    
    /* double dp2_ell_factor = 1.0 + ud.compressibility; */ /* 1.0; */
    double dp2_ell_factor = 0.5; /* 0.0; 0.5; 1.0; */
		
	const double scalefac = 1.0 / ud.Msq;
    
     const int igx = elem->igx;
     const int igy = elem->igy;
     const int igz = elem->igz;
     const int icx = elem->icx;
     const int icy = elem->icy;
     const int icz = elem->icz; 
    for(int k = igz; k < icz - igz; k++) {int l = k * icx * icy;
		for(int j = igy; j < icy - igy; j++) {int m = l + j * icx;
			for(int i = igx; i < icx - igx; i++) {int n = m + i;

                /* double p2_old    = Sol0->rhoZ[PRES][n]; */
                double p2_old       = mpv->p2_cells[n];
                double dp2_rhoY     = scalefac * (pow(Sol->rhoY[n],th.gamm) - pow(Sol0->rhoY[n],th.gamm));
                double dp2_elliptic = mpv->dp2_cells[n]; 

                Sol->rhoZ[PRES][n]  = p2_old + weight * (alpha * dp2_rhoY + (1.0-alpha) * dp2_ell_factor*dp2_elliptic);
#ifdef GRAVITY_IMPLICIT_2
                // Sol->rhoX[BUOY][n]  = Sol->rho[n] * (Sol->rho[n]/Sol->rhoY[n] - mpv->HydroState->S0[j]);
#endif
			}
		}
	}
     
    
    for (int n=0; n<elem->nc; n++) {
        /*
        double p2_old       = Sol0->rhoZ[PRES][n];
        double dp2_rhoY     = scalefac * (pow(Sol->rhoY[n],th.gamm) - pow(Sol0->rhoY[n],th.gamm));
        double dp2_elliptic = mpv->dp2_cells[n];
        
        Sol->rhoZ[PRES][n]  = (p2_old + weight * (alpha * dp2_rhoY + (1.0-alpha) * dp2_ell_factor*dp2_elliptic));
        */
        Sol->rhoZ[SOLD][n]  = Sol->rho[n]/Sol->rhoY[n];
    }
    
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double T_from_p_rho(const double p, const double rho) 
{
    extern Thermodynamic th;
    return(p/rho);
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double T_from_p_theta(const double p, const double theta) 
{
    extern Thermodynamic th;
    
    return(theta*pow(p,th.Gamma));
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double qv_sat_from_p_T(const double p, const double T) 
{
    extern Thermodynamic th;
    extern User_Data ud;
    
    double a    = th.Rg_over_Rv;
    double b    = th.Q*a/(273.16/ud.T_ref);
    double ee0  = 611/1e+05;
    double ees  = ee0*exp(b*(1-1/T));
    
    return(a * ees / (p - ees));
}


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
double compressibility(const double t) 
{    
    extern User_Data ud;
    
    switch (ud.is_compressible) {
        case 0:
            return(0.0);
            break;
        case 1:
            return(1.0);
            break;
        case -1:
        {
            /*
             double a = 0.00225;
             double b = 1.0 / 0.00225;
             return(MIN_own(1.0, MAX_own(0.0, b*(t-a))));
             */
            double a = 12.5;
            double b = 1.0 / 24.0;
            double c = MIN_own(1.0, MAX_own(0.0, b*(t-a)));
            return c;
            /*
            return(0.0);
             */
        }
            break;
        default:
            assert(0);
            break;
    }    
}




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: Eos.c,v $
 Revision 1.2  1998/03/07 09:56:43  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:30  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/

