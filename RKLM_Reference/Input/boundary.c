/*******************************************************************************
 File:   boundary.c
 Author: Nicola
 Date:   Thu Feb 19 07:13:15 CET 1998 
 *******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "Common.h"
#include "ProjectionType.h"
#include "error.h"
#include "boundary.h"
#include "userdata.h"
#include "variable.h"
#include "EOS.h"
#include "thermodynamic.h"
#include "math_own.h"
#include "memory.h"

BDRY* bdry;

static void void_min(ConsVars* Sol, const int njk, const int i);
static void void_max(ConsVars* Sol, const int njk, const int i);
static void wall_min(ConsVars* Sol, const int njk, const int i);
static void wall_max(ConsVars* Sol, const int njk, const int i);
static void slanted_wall_min(ConsVars* Sol, const int njk, const int i);
static void inflow_min(ConsVars* Sol, const int njk, const int i);
static void inflow_max(ConsVars* Sol, const int njk, const int i);
static void outflow_min(ConsVars* Sol, const int njk, const int i);
static void outflow_max(ConsVars* Sol, const int njk, const int i);
static void periodic_min(ConsVars* Sol, const int njk, const int i);
static void periodic_max(ConsVars* Sol, const int njk, const int i);
static void neumann_min(ConsVars* Sol, const int njk, const int i);
static void neumann_max(ConsVars* Sol, const int njk, const int i);
static void dirichlet_min(ConsVars* Sol, const int njk, const int i);
static void dirichlet_max(ConsVars* Sol, const int njk, const int i);
static void open_min(ConsVars* Sol, const int njk, const int i);
static void open_max(ConsVars* Sol, const int njk, const int i);

static void occupancy(
					  ConsVars* Sol, 
					  const int nijk, 
					  const int image, 
					  const int velofac);

static void (*bdry_min[])(ConsVars* sol, const int njk, const int i) = { 
void_min, 
wall_min, 
inflow_min, 
outflow_min,
periodic_min,
neumann_min,
dirichlet_min,
open_min,
slanted_wall_min
};


static void (*bdry_max[])(ConsVars* sol, const int njk, const int i) = { 
void_max, 
wall_max, 
inflow_max, 
outflow_max,
periodic_max,
neumann_max,
dirichlet_max,
open_max
};

#ifdef NODAL_PROJECTION_ONLY
static void void_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void void_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void wall_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void wall_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void slanted_wall_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void inflow_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void inflow_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void outflow_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void outflow_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void periodic_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void periodic_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void neumann_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void neumann_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void dirichlet_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void dirichlet_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void open_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);
static void open_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem);

static void (*bdry_adv_min[])(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) = { 
    void_adv_min, 
    wall_adv_min, 
    inflow_adv_min, 
    outflow_adv_min,
    periodic_adv_min,
    neumann_adv_min,
    dirichlet_adv_min,
    open_adv_min,
    slanted_wall_adv_min
};


static void (*bdry_adv_max[])(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) = { 
    void_adv_max, 
    wall_adv_max, 
    inflow_adv_max, 
    outflow_adv_max,
    periodic_adv_max,
    neumann_adv_max,
    dirichlet_adv_max,
    open_adv_max
};

static void void_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void void_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void wall_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void wall_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void slanted_wall_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void inflow_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void inflow_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void outflow_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void outflow_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void periodic_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void periodic_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void neumann_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void neumann_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void dirichlet_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void dirichlet_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void open_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void open_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);

static void (*bdry_adv_min_x[])(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) = { 
    void_adv_min_x, 
    wall_adv_min_x, 
    inflow_adv_min_x, 
    outflow_adv_min_x,
    periodic_adv_min_x,
    neumann_adv_min_x,
    dirichlet_adv_min_x,
    open_adv_min_x,
    slanted_wall_adv_min_x
};

static void (*bdry_adv_max_x[])(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) = { 
    void_adv_max_x, 
    wall_adv_max_x, 
    inflow_adv_max_x, 
    outflow_adv_max_x,
    periodic_adv_max_x,
    neumann_adv_max_x,
    dirichlet_adv_max_x,
    open_adv_max_x
};


static void void_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void void_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void wall_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void wall_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void slanted_wall_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void inflow_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void inflow_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void outflow_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void outflow_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void periodic_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void periodic_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void neumann_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void neumann_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void dirichlet_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void dirichlet_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void open_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);
static void open_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem);

static void (*bdry_adv_min_x_lat[])(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) = { 
    void_adv_min_x_lat, 
    wall_adv_min_x_lat, 
    inflow_adv_min_x_lat, 
    outflow_adv_min_x_lat,
    periodic_adv_min_x_lat,
    neumann_adv_min_x_lat,
    dirichlet_adv_min_x_lat,
    open_adv_min_x_lat,
    slanted_wall_adv_min_x_lat
};


static void (*bdry_adv_max_x_lat[])(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) = { 
    void_adv_max_x_lat, 
    wall_adv_max_x_lat, 
    inflow_adv_max_x_lat, 
    outflow_adv_max_x_lat,
    periodic_adv_max_x_lat,
    neumann_adv_max_x_lat,
    dirichlet_adv_max_x_lat,
    open_adv_max_x_lat
};


#endif

double wall_massflux(double x, double y, double wind_speed_x, double wind_speed_y);

double slanted_wall_slope(double x);
double slanted_wall_slope_3D(double x, double y);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/

void initialize_bdry(
					 const ElemSpaceDiscr* elem)
{
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int igx = elem->igx;
	const int igy = elem->igy;

    double slope_sum, slope_sum_inv;
    int i;
	
    bdry=(BDRY*)malloc(sizeof(BDRY));
	
    switch (elem->ndim) {
        case 1: {
            printf("this boundary condition setting makes no sense in 1D\n");
            exit(16);
            break;            
        }
        
        case 2: {
            bdry->wall_massflux = (double*)malloc(icx*sizeof(double));
            bdry->wall_slope    = (double*)malloc(icx*sizeof(double));
            bdry->wall_relative_slope = (double*)malloc(icx*sizeof(double));
            
            slope_sum = 0.0;
            for (i = igx; i < icx-igx; i++) {
                double slope = slanted_wall_slope(elem->x[i]);
                bdry->wall_slope[i] = slope;
                slope_sum += fabs(slope);
            }
            
            if(slope_sum <= sqrt(DBL_EPSILON)) {
                for (i = igx; i < icx-igx; i++) {
                    bdry->wall_relative_slope[i] = 0.0;
                }
            }
            else {
                slope_sum_inv = 1.0/slope_sum;
                for (i = igx; i < icx-igx; i++) {
                    bdry->wall_relative_slope[i] = slope_sum_inv * fabs(bdry->wall_slope[i]);
                }
                
                /*	*/
                slope_sum = 0.0;
                for (i = igx; i < icx-igx; i++) {
                    slope_sum += bdry->wall_relative_slope[i];
                }
                
                printf("relative slope sum = %e\n", slope_sum);
            }   
            break;
        }
        
        case 3: {
            bdry->wall_massflux = (double*)malloc(icx*icy*sizeof(double));
            bdry->wall_slope    = (double*)malloc(icx*icy*sizeof(double));
            bdry->wall_relative_slope = (double*)malloc(icx*icy*sizeof(double));
            
            slope_sum = 0.0;
            for (int j = igy; j < icy-igy; j++) {
                int nj = j*icx;
                for (int i = igx; i < icx-igx; i++) {
                    int nij = nj + i;
                    double slope = slanted_wall_slope_3D(elem->x[i], elem->y[j]);
                    bdry->wall_slope[nij] = slope;
                    slope_sum += fabs(slope);
                }
            }
            
            if(slope_sum <= sqrt(DBL_EPSILON)) {
                for (int j = igy; j < icy-igy; j++) {
                    int nj = j*icx;
                    for (int i = igx; i < icx-igx; i++) {
                        int nij = nj + i;
                        bdry->wall_relative_slope[nij] = 0.0;
                    }
                }
            }
            else {
                slope_sum_inv = 1.0/slope_sum;
                for (int j = igy; j < icy-igy; j++) {
                    int nj = j*icx;
                    for (int i = igx; i < icx-igx; i++) {
                        int nij = nj + i;
                        bdry->wall_relative_slope[nij] = slope_sum_inv * fabs(bdry->wall_slope[nij]);
                    }
                }
                
                /*	*/
                slope_sum = 0.0;
                for (int j = igy; j < icy-igy; j++) {
                    int nj = j*icx;
                    for (int i = igx; i < icx-igx; i++) {
                        int nij = nj + i;
                        slope_sum += bdry->wall_relative_slope[nij];
                    }
                }
                
                printf("relative slope sum = %e\n", slope_sum);
            }   
            break;
        } 
    
        default:
            break;
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/

void close_bdry( void )
{
    free(bdry->wall_massflux);
    free(bdry->wall_slope);
    free(bdry->wall_relative_slope);
    free(bdry);
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Bound(
		   ConsVars* Sol, 
		   const States* HydroState,
		   const double lambda, 
		   const int n, 
		   const int SplitStep) 
{
	
	/* User data */
	extern User_Data ud;
	
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	
	const int ix = elem->icx;
	const int iy = elem->icy;
	const int iz = elem->icz;
	
	int i, j, k, nk, njk, nsp;
	
	if(ud.gravity_strength[SplitStep] == 0.0) {
		/* boundary */
		for(k = 0; k < iz; k++) {
            nk = k * ix*iy;
            for(j = 0; j < iy; j++) {
                njk = nk + j * ix;
				
                /* left boundary    */
                for(i = 0 ; i < elem->igx; i++)
                    (*bdry_min[ud.bdrytype_min[SplitStep]])(Sol, njk, i);
				
				
                /* right boundaries */
                for(i = ix - elem->igx; i < ix; i++)
                    (*bdry_max[ud.bdrytype_max[SplitStep]])(Sol, njk, i);
				
            }
        }
    }
    else {
        extern Thermodynamic th;
        extern ElemSpaceDiscr* elem;
        const double g = ud.gravity_strength[SplitStep];
		const double M_LH_sq = ud.Msq;
		
        double dh = elem->dx;   /* it is "dx" because Bound() is called during OPSPLIT steps; */
        double u, v, w;
        double X[NSPEC];
	    
		/* 
		 extern double t;
		 double u_outer = velo_background(t); 
		 */
		
        for(k = 0; k < iz; k++) {
            nk = k * ix*iy;
            for(j = 0; j < iy; j++) {
                njk = nk + j * ix;
				
                /* bottom boundary ------------------------------------------------------- */				
				
                for(i=elem->igx-1;i>=0;i--) {
                    int nimage  = njk + i;
					int nlast   = nimage + 1;
                    int isource = 2*elem->igx-1 - i;
                    int nsource = njk + isource;
					
					double Z_last = Sol->rhoZ[nlast] / Sol->rho[nlast];
					double Y_last = Sol->rhoY[nlast] / Sol->rho[nlast];
					double rhou_wall;
					
                    /* copy wall-tangential velocities and scalars */
                    v = Sol->rhov[nsource] / Sol->rho[nsource]; 
                    w = Sol->rhow[nsource] / Sol->rho[nsource];
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        X[nsp] = Sol->rhoX[nsp][nsource] / Sol->rho[nsource];
                    }
                    
                    /* mirror wall-normal velocity relative to prescribed boundary mass flux */
					rhou_wall = bdry->wall_massflux[j]; 
                    u = (2.0*rhou_wall - Sol->rhou[nsource]) / Sol->rho[nsource];
					
					{
						int iimage  = i;
						/* double Yinv = (2.0*Sol->rho[nimage+1]/Sol->rhoY[nimage+1] - Sol->rho[nimage+2]/Sol->rhoY[nimage+2]); */
                        double Yinv = 1./stratification(elem->x[iimage]); 
                        double dpi  = (th.Gamma*g) * 0.5*dh*(1.0/Y_last + Yinv);
						double Z    = pow(pow(Z_last*M_LH_sq,th.Gamma) + dpi, 1.0/th.Gamma)/M_LH_sq;
                        double rhoY = pow(pow(Sol->rhoY[nlast],th.gm1) + dpi, th.gm1inv);
                        double rho  = rhoY * Yinv;
                        double p    = pow(rhoY, th.gamm);

                        /* double rhoY = HydroState->rho0[iimage]*HydroState->Y0[iimage]; */

                        Sol->rho[nimage]  = rho;
						Sol->rhou[nimage] = rho*u;
						Sol->rhov[nimage] = rho*v;
						Sol->rhow[nimage] = rho*w;
						Sol->rhoe[nimage] = rhoe(rho, u, v, w, p, Sol->geopot[nimage]);
						Sol->rhoY[nimage] = rhoY;				  /* should probably be adjusted not to take HydroState values*/
						Sol->rhoZ[nimage] = rho*Z;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Sol->rhoX[nsp][nimage] = rho*X[nsp];
                        }
					}
                }
				
                /* top boundary   ------------------------------------------------------- */
                
                for(i=0;i<elem->igx;i++) {
                    int nimage  = njk + ix - elem->igx + i;
					int nlast   = nimage - 1;
                    int isource = ix - elem->igx - 1 - i;
                    int nsource = njk + isource;
					
					double Z_last = Sol->rhoZ[nlast] / Sol->rho[nlast];
					double Y_last = Sol->rhoY[nlast] / Sol->rho[nlast];
					double rhou_wall;
					
                    /* copy wall-tangential velocities and scalars */
                    v = Sol->rhov[nsource] / Sol->rho[nsource]; 
                    w = Sol->rhow[nsource] / Sol->rho[nsource];
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        X[nsp] = Sol->rhoX[nsp][nsource] / Sol->rho[nsource];
                    }
                    
                    /* mirror wall-normal velocity */ 
					rhou_wall = 0.0;  
                    u = (2.0*rhou_wall - Sol->rhou[nsource]) / Sol->rho[nsource];
					
					{				    
						int iimage  = ix - elem->igx + i; 
                        double Yinv = 1./stratification(elem->x[iimage]); 
                        /* double Yinv = (2.0*Sol->rho[nlast]/Sol->rhoY[nlast] - Sol->rho[nlast-1]/Sol->rhoY[nlast-1]); */
                        double dpi  = -(th.Gamma*g) * 0.5*dh*(1.0/Y_last + Yinv);
                        double Z    = pow(pow(Z_last*M_LH_sq,th.Gamma) + dpi,1.0/th.Gamma)/M_LH_sq;
                        double rhoY = pow(pow(Sol->rhoY[nlast],th.gm1) + dpi, th.gm1inv);
                        double rho  = rhoY * Yinv;
                        double p    = pow(rhoY, th.gamm);

                        Sol->rho[nimage]  = rho;
						Sol->rhou[nimage] = rho*u;
						Sol->rhov[nimage] = rho*v;
						Sol->rhow[nimage] = rho*w;
						Sol->rhoe[nimage] = rhoe(rho, u, v, w, p, Sol->geopot[nimage]);
						Sol->rhoY[nimage] = rhoY;       /* should probably be adjusted not to take HydroState values*/
						Sol->rhoZ[nimage] = rho*Z;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Sol->rhoX[nsp][nimage] = rho*X[nsp];
                        }
					}
                }
            }
        }
    }
}

#ifdef NODAL_PROJECTION_ONLY
/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/

void Bound_adv_flux_x(double* rhoYu, 
                      const ElemSpaceDiscr* elem, 
                      const int SplitStep)
{
    /*
     Synchronize advective fluxes obtained from averaging from cells
     to faces with actual advective flux boundary conditions.
     */
    extern User_Data ud;
    
    switch (elem->ndim) {
        case 1:
            ERROR("Bound_adv_flux() not implemented for 1D yet.\n")
            break;
            
        case 2: 
        {
            int SplSt_lat = 1-SplitStep;
            
            (*bdry_adv_max_x[ud.bdrytype_max[SplitStep]])(rhoYu, 0, elem);
            (*bdry_adv_min_x[ud.bdrytype_min[SplitStep]])(rhoYu, 0, elem);
            
            (*bdry_adv_max_x_lat[ud.bdrytype_max[SplSt_lat]])(rhoYu, 0, elem);
            (*bdry_adv_min_x_lat[ud.bdrytype_min[SplSt_lat]])(rhoYu, 0, elem);
        }
            
            break;
            
        case 3:
            ERROR("Bound_adv_flux() not implemented for 3D yet.\n")
            break;
            
        default:
            break;
    }
}
/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Bound_adv_flux(VectorField* adv_flux, 
                    const ElemSpaceDiscr* elem)
{
    /*
     Synchronize advective fluxes obtained from averaging from cells
     to faces with actual advective flux boundary conditions.
     */
    extern User_Data ud;
    
    switch (elem->ndim) {
        case 1:
            ERROR("Bound_adv_flux() not implemented for 1D yet.\n")
            break;
            
        case 2:
            (*bdry_adv_max[ud.bdrytype_max[0]])(adv_flux, 0, elem);
            (*bdry_adv_min[ud.bdrytype_min[0]])(adv_flux, 0, elem);
            (*bdry_adv_max[ud.bdrytype_max[1]])(adv_flux, 1, elem);
            (*bdry_adv_min[ud.bdrytype_min[1]])(adv_flux, 1, elem);
            break;
            
        case 3:
            ERROR("Bound_adv_flux() not implemented for 3D yet.\n")
            break;
            
        default:
            break;
    }
    
}
#endif

#ifndef NO_BDRYCONDS_PROJ2

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void Bound_p_nodes(
                   MPV* mpv,
                   const ConsVars* Sol,
                   const ElemSpaceDiscr* elem,
                   const NodeSpaceDiscr* node,
                   const int no_of_rows) {
    
    /* User data */
    extern User_Data ud;
    
    const int icx = elem->icx;
    const int icy = elem->icy;
    
    const int inx = node->icx;
    const int iny = node->icy;
    const int inz = node->icz;

    const int igx = node->igx;
    const int igy = node->igy;
    const int igz = node->igz;
    
    const int ins[]    = {inx, iny, inz};
    const int stride[] = {1, inx, inx*iny};
    
    for (int idim = 0; idim < node->ndim; idim++) {
        
        const double g   = ud.gravity_strength[idim];

        if(g == 0) {
            /* boundary */
            for(int k = 0; k < ins[(2+idim)%3]; k++) {
                int nk = k * stride[(2+idim)%3];
                for(int j = 0; j < ins[(1+idim)%3]; j++) {
                    int njk = nk + j * stride[(1+idim)%3];
                    
                    /* left boundary */
                    for(int i = 0 ; i < no_of_rows; i++)
                        (*bdry_p_nodes_min[ud.bdrytype_min[idim]])(mpv, Sol, node, njk, i, idim);
                    
                    
                    /* right boundaries */
                    for(int i = 0; i < no_of_rows; i++)
                        (*bdry_p_nodes_max[ud.bdrytype_max[idim]])(mpv, Sol, node, njk, i, idim);
                    
                }
            }
        }
        else {
            /* rigid walls only, first row of dummy nodes only */
            
            extern Thermodynamic th;
            const double Msq = ud.Msq;
            
            double dh = node->dxyz[idim];  /* it is "dxyz[idim]" because Bound_p_nodes() is called outside OPSPLIT steps; */
            
            for(int k = igz; k < inz-igz; k++) {
                int ln = k * inx*iny;
                int lc = k * icx*icy;
                for(int i = igx; i < inx-igx; i++) {
                    int mn = ln + i;
                    int mc = lc + i;
                    
                    /* bottom boundary ------------------------------------------------------- */
                    
                    for(int j = igy-1; j > 0; j--) {
                        int nnimage  = mn + j*inx;
                        int nnlast   = nnimage + inx;
                        int jnsource = igy + (igy - j);
                        int nnsource = mn + jnsource*inx;
                        int ncp      = mc + j*icx;
                        int ncm      = mc + j*icx - 1;
                        int ncp_last = mc + (jnsource-1)*icx;
                        int ncm_last = mc + (jnsource-1)*icx - 1;

                        double p_last    = mpv->p2_nodes[nnlast];
                        double p_source  = mpv->p2_nodes[nnsource];
                        double dp_source = mpv->dp2_nodes[nnsource];
                        double Yinv      = 0.5*(Sol->rho[ncm] / Sol->rhoY[ncm] + Sol->rho[ncp] / Sol->rhoY[ncp]);
                        double Yinv_last = 0.5*(Sol->rho[ncm_last] / Sol->rhoY[ncm_last] + Sol->rho[ncp_last] / Sol->rhoY[ncp_last]);
                        double phy_image  = pow(pow(p_last*Msq,th.Gamma) + (th.Gamma*g) * dh * Yinv, th.Gammainv)/Msq;
                        double phy_source = pow(pow(p_last*Msq,th.Gamma) - (th.Gamma*g) * dh * Yinv_last, th.Gammainv)/Msq;

                        mpv->p2_nodes[nnimage]  = phy_image + (p_source-phy_source);
                        mpv->dp2_nodes[nnimage] = dp_source;
                    }
                    
                    /* top boundary   ------------------------------------------------------- */
                    
                    for(int j = iny-igy; j < iny-1; j++) {
                        int nnimage  = mn + j*inx;
                        int nnlast   = nnimage - inx;
                        int jnsource = (iny-igy-1) - (j - (iny-igy-1));
                        int nnsource = mn + jnsource*inx;
                        int ncp      = mc + j*icx;
                        int ncm      = mc + j*icx - 1;
                        int ncp_last = mc + jnsource*icx;
                        int ncm_last = mc + jnsource*icx - 1;

                        double p_last     = mpv->p2_nodes[nnlast];
                        double p_source   = mpv->p2_nodes[nnsource];
                        double dp_source  = mpv->dp2_nodes[nnsource];
                        double Yinv       = 0.5*(Sol->rho[ncm] / Sol->rhoY[ncm] + Sol->rho[ncp] / Sol->rhoY[ncp]);
                        double Yinv_last  = 0.5*(Sol->rho[ncm_last] / Sol->rhoY[ncm_last] + Sol->rho[ncp_last] / Sol->rhoY[ncp_last]);
                        double phy_image  = pow(pow(p_last*Msq,th.Gamma) - (th.Gamma*g) * dh * Yinv, th.Gammainv)/Msq;
                        double phy_source = pow(pow(p_last*Msq,th.Gamma) + (th.Gamma*g) * dh * Yinv_last, th.Gammainv)/Msq;
                        
                        mpv->p2_nodes[nnimage]  = phy_image + (p_source-phy_source);
                        mpv->dp2_nodes[nnimage] = dp_source;
                    }
                }
            }
        }
    }
}

#endif /* NO_BDRYCONDS_PROJ2 */

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void set_wall_massflux(
					   BDRY* bdry, 
					   const ConsVars* Sol0, 
					   const ElemSpaceDiscr* elem)
{
	const int icx = elem->icx;
	const int icy = elem->icy;
	const int igx = elem->igx;
	const int igy = elem->igy;
	double wall_mf, wall_flux_balance;
	int i;
	
    switch (elem->ndim) {
        case 1: {
            printf("\nwall flux in 1D makes no sense\n");
            exit(15);
            break;                
        }
            
        case 2: {
            /* first guess for wall mass fluxes */
            const int nstart = elem->igy*elem->icx;

            wall_flux_balance = 0.0;
            for (i = igx; i < icx-igx; i++) {
                int n = nstart + i;
                wall_mf = wall_massflux(elem->x[i], elem->y[i], Sol0->rhou[n]/Sol0->rho[n], Sol0->rhov[n]/Sol0->rho[n]);
                bdry->wall_massflux[i] = wall_mf;
                wall_flux_balance += wall_mf;
            }
            
            for(i=0; i<elem->igx; i++) {
                bdry->wall_massflux[i] = bdry->wall_massflux[icx-1-i] = 0.0; 
            }
            
            /* correction for zero net flux */
            for (i = igx; i < icx-igx; i++) {
                bdry->wall_massflux[i] -= wall_flux_balance * bdry->wall_relative_slope[i];
            }
            
            /*	*/
            {
                double flux_sum = 0.0;
                for (i = elem->igx; i < elem->icx-elem->igx; i++) {
                    flux_sum += bdry->wall_massflux[i];
                }
                
                printf("wall flux sum = %e\n", flux_sum);
            }
            break;
        }
            
        case 3: {
            const int nstart = elem->igz*elem->icx*elem->icy;

            /* first guess for wall mass fluxes before total flux correction (for elliptic solvability cond.) */
            wall_flux_balance = 0.0;
            for (int j = igy; j < icy-igy; j++) {
                int njk = nstart + j*icx;
                int nj  = j*icx;
                for (int i = igx; i < icx-igx; i++) {
                    int nijk = njk + i;
                    int nij  = nj  + i;
                    wall_mf = wall_massflux(elem->x[i], elem->y[i], Sol0->rhou[nijk]/Sol0->rho[nijk], Sol0->rhov[nijk]/Sol0->rho[nijk]);
                    bdry->wall_massflux[nij] = wall_mf;
                    wall_flux_balance += wall_mf;
                }
            }
            
            /* set wall flux in dummy cells to zero */
            for(int j=0; j<icy; j++) {
                int nj = j * icx;
                for(int i=0; i<igx; i++) {
                    int nij_left  = nj + i;
                    int nij_right = nj + icx-1-i;
                    bdry->wall_massflux[nij_left] = bdry->wall_massflux[nij_right] = 0.0; 
                }
            }            

            for(int i=0; i<icx; i++) {
                int ni = i;
                for(int j=0; j<igy; j++) {
                    int nij_left  = ni + j*icx;
                    int nij_right = ni + (icy-1-j)*icx;
                    bdry->wall_massflux[nij_left] = bdry->wall_massflux[nij_right] = 0.0; 
                }
            }            
            
            /* correction for zero net flux (needed for elliptic solvability) */
            for (int j = igy; j < icy-igy; j++) {
                int nj = j * icx;
                for (int i = igx; i < icx-igx; i++) {
                    int nij  = nj + i;
                    bdry->wall_massflux[nij] -= wall_flux_balance * bdry->wall_relative_slope[nij];
                }
            }
            
            /*	*/
            {
                double flux_sum = 0.0;
                for (int j = igy; j < icy-igy; j++) {
                    int nj = j * icx;
                    for (int i = igx; i < icx-igx; i++) {
                        int nij  = nj + i;
                        flux_sum += bdry->wall_massflux[nij];
                    }
                }
                
                printf("wall flux sum = %e\n", flux_sum);
            }
            
            break;
        }
            
    }
}

static void void_min(ConsVars* Sol, const int njk, const int i) {}

static void void_max(ConsVars* Sol, const int njk, const int i) {}

static void wall_min(ConsVars* Sol, const int njk, const int i) 
{
	extern ElemSpaceDiscr* elem;
	int velofac; /* Factor for sign of velocity */
	int image, nijk;
	
	nijk    = njk + elem->igx-1 - i; 
	image   = njk + elem->igx + i; 
	
	velofac = -1.0;
	
	occupancy( Sol, nijk, image, velofac );
}

static void slanted_wall_min(ConsVars* Sol, const int njk, const int i) 
{
	/* 
     Very specialized implementation of perturbational curved
     wall boundary condition.
	 
     Goal is to run a test case for anelastic flow over a shallow
	 hill, which should induce gravity waves. 
	 
	 I only take into account the "y-direction", and only the "min"-
     boundary. 
	 */
	
	extern ElemSpaceDiscr* elem;
	const int ix = elem->icx;
	const int iy = elem->icy;
	
	double rhou_wall, dhdx;
	int image, nijk, jj, kk;
	
	nijk    = njk + elem->igx-1 - i; 
	image   = njk + elem->igx + i; 
	
	/* extract horizontal position of current y-column */
	kk = nijk / (iy*ix);
	jj = (nijk - iy*ix*kk) / ix;
	
	/* provide wall slope */
	dhdx = slanted_wall_slope(elem->y[jj]);
	
	rhou_wall = dhdx * Sol->rhov[image];
	
	Sol->rho[nijk]    = Sol->rho[image];
	Sol->rhou[nijk]   = rhou_wall - (Sol->rhou[image]-rhou_wall);
	Sol->rhov[nijk]   = Sol->rhov[image];
	Sol->rhow[nijk]   = Sol->rhow[image];
	Sol->rhoe[nijk]   = Sol->rhoe[image];
	Sol->rhoY[nijk]   = Sol->rhoY[image];
	Sol->rhoZ[nijk]   = Sol->rhoZ[image];
}

static void wall_max(ConsVars* Sol, const int njk, const int i)
{  
	
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	const int ix = elem->icx;
	
	int velofac; /* Factor for sign of velocity */
	int image, nijk;
	
	nijk    = njk + i;
	image   = njk + (ix-elem->igx-1) - (i-(ix-elem->igx));
	
	velofac = -1.0;
	
	occupancy( Sol, nijk, image, velofac );
}

static void inflow_min(ConsVars* Sol, const int njk, const int i) {
	ERROR("function not available");
}

static void inflow_max(ConsVars* Sol, const int njk, const int i) {
	ERROR("function not available");
}

static void outflow_min(ConsVars* Sol, const int njk, const int i) {
	ERROR("function not available");
}

static void outflow_max(ConsVars* Sol, const int njk, const int i) {
	ERROR("function not available");
}

static void periodic_min(ConsVars* Sol, const int njk, const int i)
{
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	const int ix = elem->icx;
	
	int velofac; /* Factor for sign of velocity */
	int image, nijk;
	
	nijk    = njk + i;
	image   = nijk + (ix-2*elem->igx);
	
	velofac = 1.0;
	
	occupancy( Sol, nijk, image, velofac );
}

static void periodic_max(ConsVars* Sol, const int njk, const int i)
{
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	const int ix = elem->icx;
	
	int velofac; /* Factor for sign of velocity */
	int image, nijk;
	
	nijk    = njk + i;
	image   = nijk - (ix-2*elem->igx);
	
	velofac = 1.0;
	
	occupancy( Sol, nijk, image, velofac );
}

static void neumann_min(ConsVars* Sol, const int njk, const int i)
{
	extern ElemSpaceDiscr* elem;
	int velofac; /* Factor for sign of velocity */ 
	int image, nijk;
	
	nijk    = njk + elem->igx-1 - i;
	image   = njk + elem->igx + i;
	
	velofac = 1.0;
	
	occupancy( Sol, nijk, image, velofac );
}

static void neumann_max(ConsVars* Sol, const int njk, const int i)
{
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	const int ix = elem->icx;
	
	int velofac; /* Factor for sign of velocity */
	int image, nijk;
	
	nijk    = njk + i;
	image   = njk + (ix-elem->igx-1) - (i-(ix-elem->igx));
	
	velofac = 1.0;
	
	occupancy( Sol, nijk, image, velofac );
}

static void dirichlet_min(ConsVars* Sol, const int njk, const int i)
{
	extern MPV* mpv;
	extern ElemSpaceDiscr* elem;
	extern double t;
	
	int nijk;
	
	int k = njk / (elem->icy*elem->icx);
	int j = (njk - k*elem->icy*elem->icx) / elem->icx;
	
	double p_outer   = mpv->HydroState->p0[j];
	double rho_outer = mpv->HydroState->rho0[j];
	double S2_outer  = mpv->HydroState->Y0[j];
	double p2_outer  = mpv->HydroState->p20[j];
	double u_outer   = velo_background(t);
	
	nijk    = njk + elem->igx-1 - i;
	
	Sol->rho[nijk]    = rho_outer;
	Sol->rhou[nijk]   = rho_outer * u_outer;
	Sol->rhov[nijk]   = 0.0;
	Sol->rhow[nijk]   = 0.0;
	Sol->rhoe[nijk]   = rhoe(rho_outer, u_outer, 0.0, 0.0, p_outer, Sol->geopot[nijk]);
	Sol->rhoY[nijk]   = rho_outer * S2_outer;
	Sol->rhoZ[nijk]   = rho_outer * p2_outer;  
}

static void dirichlet_max(ConsVars* Sol, const int njk, const int i)
{
	extern MPV* mpv;
	extern ElemSpaceDiscr* elem;
	extern double t;
	
	int nijk;
	
	int k = njk / (elem->icy*elem->icx);
	int j = (njk - k*elem->icy*elem->icx) / elem->icx;
	
	double p_outer   = mpv->HydroState->p0[j];
	double rho_outer = mpv->HydroState->rho0[j];
	double S2_outer  = mpv->HydroState->Y0[j];
	double p2_outer  = mpv->HydroState->p20[j];
	double u_outer   = velo_background(t);
	
	nijk = njk + i;
	
	Sol->rho[nijk]    = rho_outer;
	Sol->rhou[nijk]   = rho_outer * u_outer;
	Sol->rhov[nijk]   = 0.0;
	Sol->rhow[nijk]   = 0.0;
	Sol->rhoe[nijk]   = rhoe(rho_outer, u_outer, 0.0, 0.0, p_outer, Sol->geopot[nijk]);
	Sol->rhoY[nijk]   = rho_outer * S2_outer;
	Sol->rhoZ[nijk]   = rho_outer * p2_outer;
}

static void open_min(ConsVars* Sol, const int njk, const int i)
{
	extern ElemSpaceDiscr* elem;
	int velofac; /* Factor for sign of velocity */ 
	int image, nijk;
	
	nijk    = njk + elem->igx-1 - i;
	image   = njk + elem->igx + i;
	
	velofac = 1.0;
	
	occupancy( Sol, nijk, image, velofac );
}

static void open_max(ConsVars* Sol, const int njk, const int i)
{
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	const int ix = elem->icx;
	
	int velofac; /* Factor for sign of velocity */
	int image, nijk;
	
	nijk    = njk + i;
	image   = njk + (ix-elem->igx-1) - (i-(ix-elem->igx));
	
	velofac = 1.0;
	
	occupancy( Sol, nijk, image, velofac );
}

#ifdef NODAL_PROJECTION_ONLY
static void void_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {}
static void void_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {}

static void wall_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) 
{
    const int icx = elem->icx;
    const int icy = elem->icy;

    const int ifx = elem->ifx;
    const int ify = elem->ify;

    const int igx = elem->igx;
    const int igy = elem->igy;

    if (idim >= elem->ndim) {
        return;
    }

    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }

    switch (idim) {
        case 0:
            for (int j=0; j<icy; j++) {
                int nj   = j*ifx;
                int njim = nj+igx;
                adv_flux->x[njim] = 0.0;
            }
            break;
            
        case 1:
            for (int i=0; i<icx; i++) {
                int ni   = i*ify;
                int nijm = ni+igy;
                adv_flux->y[nijm] = 0.0;
            }
            break;
                        
        default:
            break;
    }
}

static void slanted_wall_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void wall_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem)
{  
    const int icx = elem->icx;
    const int icy = elem->icy;
    
    const int ifx = elem->ifx;
    const int ify = elem->ify;
    
    const int igx = elem->igx;
    const int igy = elem->igy;
    
    if (idim >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    switch (idim) {
        case 0:
            for (int j=0; j<icy; j++) {
                int nj   = j*ifx;
                int njip = nj+ifx-igx-1;
                adv_flux->x[njip] = 0.0;
            }
            break;
            
        case 1:
            for (int i=0; i<icx; i++) {
                int ni   = i*ify;
                int nijp = ni+ify-igy-1;
                adv_flux->y[nijp] = 0.0;
            }
            break;
            
        default:
            break;
    }
}

static void inflow_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void inflow_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void outflow_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void outflow_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void periodic_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem)
{
    const int icx = elem->icx;
    const int icy = elem->icy;
    
    const int ifx = elem->ifx;
    const int ify = elem->ify;
    
    const int igx = elem->igx;
    const int igy = elem->igy;
    
    if (idim >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    switch (idim) {
        case 0:
            for (int j=0; j<icy; j++) {
                int nj   = j*ifx;
                int njim = nj+igx;
                int njip = nj+ifx-igx-1;
                adv_flux->x[njim] = adv_flux->x[njip];
            }
            break;
            
        case 1:
            for (int i=0; i<icx; i++) {
                int ni   = i*ify;
                int nijm = ni+igy;
                int nijp = ni+ify-igy-1;
                adv_flux->y[nijm] = adv_flux->y[nijp];
            }
            break;
            
        default:
            break;
    }    
}

static void periodic_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem)
{
    const int icx = elem->icx;
    const int icy = elem->icy;
    
    const int ifx = elem->ifx;
    const int ify = elem->ify;
    
    const int igx = elem->igx;
    const int igy = elem->igy;
    
    if (idim >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    switch (idim) {
        case 0:
            for (int j=0; j<icy; j++) {
                int nj   = j*ifx;
                int njim = nj+igx;
                int njip = nj+ifx-igx-1;
                adv_flux->x[njip] = adv_flux->x[njim];
            }
            break;
            
        case 1:
            for (int i=0; i<icx; i++) {
                int ni   = i*ify;
                int nijm = ni+igy;
                int nijp = ni+ify-igy-1;
                adv_flux->y[nijp] = adv_flux->y[nijm];
            }
            break;
            
        default:
            break;
    }    
}

static void neumann_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void neumann_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void dirichlet_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void dirichlet_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void open_adv_min(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void open_adv_max(VectorField* adv_flux, const int idim, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

/* ======================================================================== */

static void void_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {}
static void void_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {}

static void wall_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) 
{
    const int icy = elem->icy;    
    const int ifx = elem->ifx;    
    const int igx = elem->igx;
    
    if (SplitStep >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    for (int j=0; j<icy; j++) {
        int nj   = j*ifx;
        int njim = nj+igx;
        adv_flux_x[njim] = 0.0;
    }
}

static void slanted_wall_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void wall_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem)
{  
    const int icy = elem->icy;    
    const int ifx = elem->ifx;    
    const int igx = elem->igx;
    
    if (SplitStep >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    for (int j=0; j<icy; j++) {
        int nj   = j*ifx;
        int njip = nj+ifx-igx-1;
        adv_flux_x[njip] = 0.0;
    }
}

static void inflow_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void inflow_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void outflow_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void outflow_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void periodic_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem)
{
    const int icy = elem->icy;
    const int ifx = elem->ifx;
    const int igx = elem->igx;
    
    if (SplitStep >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    for (int j=0; j<icy; j++) {
        int nj   = j*ifx;
        int njim = nj+igx;
        int njip = nj+ifx-igx-1;
        adv_flux_x[njim] = adv_flux_x[njip];
    }
}

static void periodic_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem)
{
    const int icy = elem->icy;    
    const int ifx = elem->ifx;    
    const int igx = elem->igx;
    
    if (SplitStep >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    for (int j=0; j<icy; j++) {
        int nj   = j*ifx;
        int njim = nj+igx;
        int njip = nj+ifx-igx-1;
        adv_flux_x[njip] = adv_flux_x[njim];
    }
}

static void neumann_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void neumann_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void dirichlet_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void dirichlet_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void open_adv_min_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void open_adv_max_x(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

/* ======================================================================== */

static void void_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {}
static void void_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {}

static void wall_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) 
{
    /* 
     lateral-to-x boundary is a rigid wall 
     */
    const int ifx = elem->ifx;    
    const int igy = elem->igy;
    
    if (SplitStep >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    for (int j=0; j<igy; j++) {
        int nimage  = (igy-1-j)*ifx;
        int nsource = (igy+j)*ifx;
        for (int i=0; i<ifx; i++) {
            int mimage  = nimage  + i;
            int msource = nsource + i;
            adv_flux_x[mimage] = adv_flux_x[msource];
        }
    }
}

static void slanted_wall_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void wall_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem)
{  
    const int icy = elem->icy;    
    const int ifx = elem->ifx;    
    const int igy = elem->igy;
    
    if (SplitStep >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    for (int j=0; j<igy; j++) {
        int nimage  = (icy - igy + j)*ifx;
        int nsource = (icy - igy - 1 - j)*ifx;
        for (int i=0; i<ifx; i++) {
            int mimage  = nimage  + i;
            int msource = nsource + i;
            adv_flux_x[mimage] = adv_flux_x[msource];
        }
    }
}

static void inflow_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void inflow_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void outflow_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void outflow_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}

static void periodic_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem)
{
    const int icy = elem->icy;
    const int ifx = elem->ifx;
    const int igy = elem->igy;
    
    if (SplitStep >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    for (int j=0; j<igy; j++) {
        int nimage  = (igy-1-j)*ifx;
        int nsource = (icy-igy-1-j)*ifx;
        for (int i=0; i<ifx; i++) {
            int mimage  = nimage  + i;
            int msource = nsource + i;
            adv_flux_x[mimage] = adv_flux_x[msource];
        }
    }
}

static void periodic_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem)
{
    const int icy = elem->icy;    
    const int ifx = elem->ifx;    
    const int igy = elem->igy;
    
    if (SplitStep >= elem->ndim) {
        return;
    }
    
    if (elem->ndim != 2) {
        ERROR("adv_flux boundary conditions implemented in 2D only so far.");
    }
    
    for (int j=0; j<igy; j++) {
        int nimage  = (icy-igy+j)*ifx;
        int nsource = (igy+j)*ifx;
        for (int i=0; i<ifx; i++) {
            int mimage  = nimage  + i;
            int msource = nsource + i;
            adv_flux_x[mimage] = adv_flux_x[msource];
        }
    }
}

static void neumann_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void neumann_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void dirichlet_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void dirichlet_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void open_adv_min_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


static void open_adv_max_x_lat(double* adv_flux_x, const int SplitStep, const ElemSpaceDiscr* elem) {
    ERROR("function not available");
}


#endif

/* ======================================================================== */
/* ======================================================================== */


static void occupancy(
					  ConsVars* Sol, 
					  const int nijk, 
					  const int image, 
					  const int velofac) {
    extern User_Data ud;
    int nsp;
	
	Sol->rho[nijk]    = Sol->rho[image];
	Sol->rhou[nijk]   = Sol->rhou[image] * velofac;
	Sol->rhov[nijk]   = Sol->rhov[image];
	Sol->rhow[nijk]   = Sol->rhow[image];
	Sol->rhoe[nijk]   = Sol->rhoe[image];
	Sol->rhoY[nijk]   = Sol->rhoY[image];
	Sol->rhoZ[nijk]   = Sol->rhoZ[image];
    for (nsp = 0; nsp < ud.nspec; nsp++) {
        Sol->rhoX[nsp][nijk]   = Sol->rhoX[nsp][image];
    }
}

double slanted_wall_slope(double x){
	extern User_Data ud;
	double hill_height  = ud.hill_height;
	double length_scale_inv = 1.0/ud.hill_length_scale; 
	double x_sc = x*length_scale_inv;
	
	return(- hill_height * 2.0*x_sc / ((1.0 + x_sc*x_sc)*(1.0 + x_sc*x_sc)) * length_scale_inv);
	
}

double slanted_wall_slope_3D(double x, double y){
	extern User_Data ud;
	double hill_height  = ud.hill_height;
	double length_scale_inv = 1.0/ud.hill_length_scale; 
	double x_sc = x*length_scale_inv;
	
	return(- hill_height * 2.0*x_sc / ((1.0 + x_sc*x_sc)*(1.0 + x_sc*x_sc)) * length_scale_inv);
	
}

double velo_background(double t){
	
	extern User_Data ud;
	
	return(ud.wind_speed);
	
	/*
	 double t_1 = 0.1;
	 double u_max_infty = 0.00;
	 
	 return((t < t_1 ? u_max_infty * 0.5 * (1.0 - cos(PI*t/t_1)) : u_max_infty));
	 */
}

double wall_massflux(double x, double y, double wind_speed_x, double wind_speed_y){
	extern User_Data ud;
    
	double hill_height  = ud.hill_height;
	double length_scale_inv = 1.0/ud.hill_length_scale; 
	double x_sc = x*length_scale_inv;
	
    /*
	return(- wind_speed * hill_height * 2.0*x_sc / ((1.0 + x_sc*x_sc)*(1.0 + x_sc*x_sc)) * length_scale_inv);
     */
    return(- wind_speed_x * hill_height * 2.0*x_sc / ((1.0 + x_sc*x_sc)*(1.0 + x_sc*x_sc)) * length_scale_inv
           - wind_speed_y * 0.0 );
}

/* ========================================================================= */

void check_flux_bcs(
					States* Lefts, 
					States* Rights,
					const int nmax,
					const int kcache,
					const int njump,
					ElemSpaceDiscr* elem,
					const int SplitStep) {
	
	extern User_Data ud;
	
	const double g = ud.gravity_strength[ud.gravity_direction];
	
    int i, nfull, ifull, jfull, nsp; 
    
	switch (SplitStep) {
		case 1: {
			if(ud.bdrytype_min[SplitStep] == WALL) {
				for(i=elem->igx; i<nmax-elem->igx+1; i++){
					nfull = kcache*njump + i;
					jfull = nfull/elem->icx;
					ifull = nfull%elem->icx;
					if(ifull==elem->igx) {
#ifdef PERTURBED_WALL
						double rhou_wall = bdry->wall_massflux[jfull];
						Lefts->rhou[i-1] = Rights->rhou[i] = rhou_wall;
						Lefts->rhoY[i-1] = Rights->rhoY[i] = Rights->rho[i] * stratification(0.0);
#else
						Lefts->rhou[i-1] = -Rights->rhou[i];
						Lefts->rhoY[i-1] = Rights->rhoY[i];
#endif
						
						Lefts->rho[i-1]  = Rights->rho[i];
						Lefts->rhov[i-1] = Rights->rhov[i];
						Lefts->rhow[i-1] = Rights->rhow[i];
						Lefts->rhoe[i-1] = Rights->rhoe[i];
						Lefts->rhoZ[i-1] = Rights->rhoZ[i];
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Lefts->rhoX[nsp][i-1] = Rights->rhoX[nsp][i];
                        }
					}
				}
			}
			
			if(ud.bdrytype_max[SplitStep] == WALL) {
				for(i=elem->igx; i<nmax - elem->igx+1; i++){
					nfull = kcache*njump + i;
					ifull = nfull%elem->icx;
					if(ifull == elem->icx - elem->igx) {
						Rights->rho[i] = Lefts->rho[i-1];
						Rights->rhou[i] = - Lefts->rhou[i-1];
						Rights->rhov[i] = Lefts->rhov[i-1];
						Rights->rhow[i] = Lefts->rhow[i-1];
						Rights->rhoe[i] = Lefts->rhoe[i-1];
						Rights->rhoY[i] = Lefts->rhoY[i-1];
						Rights->rhoZ[i] = Lefts->rhoZ[i-1];
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Rights->rhoX[nsp][i] = Lefts->rhoX[nsp][i-1];
                        }
					}
				}
			}
			break;
		}
		default: {
			if(ud.bdrytype_min[SplitStep] == WALL) {
				for(i=elem->igx; i<nmax-elem->igx; i++){
					nfull = kcache*njump + i;
					ifull = nfull%elem->icx;
					if(ifull==elem->igx) {
						Lefts->rho[i-1]  = Rights->rho[i];
						Lefts->rhou[i-1] = - Rights->rhou[i];
						Lefts->rhov[i-1] = Rights->rhov[i];
						Lefts->rhow[i-1] = Rights->rhow[i];
						Lefts->rhoe[i-1] = Rights->rhoe[i];
						Lefts->rhoY[i-1] = Rights->rhoY[i];
						Lefts->rhoZ[i-1] = Rights->rhoZ[i];
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Lefts->rhoX[nsp][i-1] = Rights->rhoX[nsp][i];
                        }
					}
				}
			}
			else if(ud.bdrytype_min[SplitStep] == DIRICHLET) {
				extern MPV* mpv;
				extern User_Data ud;
				for(i=elem->igx; i<nmax-elem->igx; i++){
					nfull = kcache*njump + i;
					jfull = nfull/elem->icx;
					ifull = nfull%elem->icx;
					if(ifull==elem->igx) {
						double rho = mpv->HydroState->rho0[jfull];
						double p   = mpv->HydroState->p0[jfull];
						double y   = elem->y[jfull];
						double S2  = mpv->HydroState->Y0[jfull];
						double p2  = mpv->HydroState->p20[jfull];
						
						Rights->rho[i]  = Lefts->rho[i-1]  = rho;
						Rights->rhou[i] = Lefts->rhou[i-1] = rho * ud.wind_speed;
						Rights->rhov[i] = Lefts->rhov[i-1] = 0.0;
						Rights->rhow[i] = Lefts->rhow[i-1] = 0.0;
						Rights->rhoe[i] = Lefts->rhoe[i-1] = rhoe(rho, ud.wind_speed, 0.0, 0.0, p, g*y);
						Rights->rhoY[i] = Lefts->rhoY[i-1] = rho * S2; 
						Rights->rhoZ[i] = Lefts->rhoZ[i-1] = rho * p2;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Rights->rhoX[nsp][i] = Lefts->rhoX[nsp][i-1] = 99999.9;
                        }
					}
				}
			}
			
			if(ud.bdrytype_max[SplitStep] == WALL) {
				for(i=elem->igx; i<nmax - elem->igx; i++){
					nfull = kcache*njump + i;
					ifull = nfull%elem->icx;
					if(ifull == elem->icx - elem->igx) {
						Rights->rho[i] = Lefts->rho[i-1];
						Rights->rhou[i] = - Lefts->rhou[i-1];
						Rights->rhov[i] = Lefts->rhov[i-1];
						Rights->rhow[i] = Lefts->rhow[i-1];
						Rights->rhoe[i] = Lefts->rhoe[i-1];
						Rights->rhoY[i] = Lefts->rhoY[i-1];
						Rights->rhoZ[i] = Lefts->rhoZ[i-1];
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Rights->rhoX[nsp][i] = Lefts->rhoX[nsp][i-1];
                        }
					}
				}
			}
			else if(ud.bdrytype_max[SplitStep] == DIRICHLET) {
				extern MPV* mpv;
				for(i=elem->igx; i<nmax - elem->igx; i++){
					nfull = kcache*njump + i;
					jfull = nfull/elem->icx;
					ifull = nfull%elem->icx;
					if(ifull == elem->icx - elem->igx) {
						double rho = mpv->HydroState->rho0[jfull];
						double p   = mpv->HydroState->p0[jfull];
						double y   = elem->y[jfull];
						double S2  = mpv->HydroState->Y0[jfull];
						double p2  = mpv->HydroState->p20[jfull];
						
						Lefts->rho[i-1]  = Rights->rho[i]  = rho;
						Lefts->rhou[i-1] = Rights->rhou[i] = rho * ud.wind_speed;
						Lefts->rhov[i-1] = Rights->rhov[i] = 0.0;
						Lefts->rhow[i-1] = Rights->rhow[i] = 0.0;
						Lefts->rhoe[i-1] = Rights->rhoe[i] = rhoe(rho, ud.wind_speed, 0.0, 0.0, p, g*y);
						Lefts->rhoY[i-1] = Rights->rhoY[i] = rho * S2;
						Lefts->rhoZ[i-1] = Rights->rhoZ[i] = rho * p2;
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Lefts->rhoX[nsp][i-1] = Rights->rhoX[nsp][i] = 9999.9;
                        }
					}
				}
			}
			break;
		}
	}
}

/* ============================================================================= */

static void (*rotate[])(ConsVars* Sol, double* rhs, double *Yinvbg, const enum Direction dir) = {NULL, rotate2D, rotate3D};

void Set_Explicit_Boundary_Data(
                                ConsVars* Sol,
                                const ElemSpaceDiscr* elem,
                                const MPV* mpv) 
{
    extern double *Yinvbg;
    int SplitStep;
    
    for(SplitStep = 0; SplitStep < elem->ndim; SplitStep++) { 
        const double lambda = 1.0;
        Bound(Sol, mpv->HydroState, lambda, elem->nc, SplitStep); 
        if(SplitStep < elem->ndim - 1) (*rotate[elem->ndim - 1])(Sol, mpv->Level[0]->rhs, Yinvbg, FORWARD);
    }         
    /* rotate back */          
    for(SplitStep = elem->ndim-1; SplitStep > 0; SplitStep--) {
        (*rotate[elem->ndim - 1])(Sol, mpv->Level[0]->rhs, Yinvbg, BACKWARD);
    }
    
}


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: boundary.c,v $
 Revision 1.1  1998/03/01 18:43:32  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
