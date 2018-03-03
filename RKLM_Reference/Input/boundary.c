/*******************************************************************************
 File:   boundary.c
 Author: Nicola
 Date:   Thu Feb 19 07:13:15 CET 1998 
 *******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "Common.h"
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
		   const double lambda, 
		   const int n, 
		   const int SplitStep, 
           const int setZ) 
{
	
	/* User data */
	extern User_Data ud;
    extern MPV* mpv;

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
        const int compressible = ud.is_compressible;
		
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
						/* double S = (2.0*Sol->rho[nimage+1]/Sol->rhoY[nimage+1] - Sol->rho[nimage+2]/Sol->rhoY[nimage+2]); */
                        double S    = 1./stratification(elem->x[iimage]); 
                        double dpi  = (th.Gamma*g) * 0.5*dh*(1.0/Y_last + S);
                        double rhoY = (compressible == 1 ? pow(pow(Sol->rhoY[nlast],th.gm1) + dpi, th.gm1inv) : mpv->HydroState->rhoY0[i]);
                        double rho  = rhoY * S;
                        double p    = pow(rhoY, th.gamm);

                        /* treat as p/Pbar
                         not implemented - but worth a try
                         */

                        Sol->rho[nimage]  = rho;
						Sol->rhou[nimage] = rho*u;
						Sol->rhov[nimage] = rho*v;
						Sol->rhow[nimage] = rho*w;
						Sol->rhoe[nimage] = rhoe(rho, u, v, w, p, Sol->geopot[nimage]);
						Sol->rhoY[nimage] = rhoY;				  /* should probably be adjusted not to take HydroState values*/
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
                        double S = 1./stratification(elem->x[iimage]); 
                        /* double S = (2.0*Sol->rho[nlast]/Sol->rhoY[nlast] - Sol->rho[nlast-1]/Sol->rhoY[nlast-1]); */
                        double dpi  = -(th.Gamma*g) * 0.5*dh*(1.0/Y_last + S);
                        double rhoY = (compressible == 1 ? pow(pow(Sol->rhoY[nlast],th.gm1) + dpi, th.gm1inv) : mpv->HydroState->rhoY0[iimage]);
                        double rho  = rhoY * S;
                        double p    = pow(rhoY, th.gamm);

                        /* treat as  p / Pbar 
                         not implemented - but worth trying 
                         */

                        Sol->rho[nimage]  = rho;
						Sol->rhou[nimage] = rho*u;
						Sol->rhov[nimage] = rho*v;
						Sol->rhow[nimage] = rho*w;
						Sol->rhoe[nimage] = rhoe(rho, u, v, w, p, Sol->geopot[nimage]);
						Sol->rhoY[nimage] = rhoY;       /* should probably be adjusted not to take HydroState values*/
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            Sol->rhoX[nsp][nimage] = rho*X[nsp];
                        }
					}
                }
            }
        }
    }
}

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
                        double rhou_wall = bdry->wall_massflux[jfull];
						Lefts->rhou[i-1] = Rights->rhou[i] = rhou_wall;
						Lefts->rhoY[i-1] = Rights->rhoY[i] = Rights->rho[i] * stratification(0.0);
						
						Lefts->rho[i-1]  = Rights->rho[i];
						Lefts->rhov[i-1] = Rights->rhov[i];
						Lefts->rhow[i-1] = Rights->rhow[i];
						Lefts->rhoe[i-1] = Rights->rhoe[i];
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

static void (*rotate[])(ConsVars* Sol, double* rhs, double *Sbg, double *buoyS, const enum Direction dir) = {NULL, rotate2D, rotate3D};

void Set_Explicit_Boundary_Data(
                                ConsVars* Sol,
                                const ElemSpaceDiscr* elem,
                                const MPV* mpv,
                                const int setZ) 
{
    extern double *Sbg;
    extern double *buoyS;
    int SplitStep;
    
    for(SplitStep = 0; SplitStep < elem->ndim; SplitStep++) { 
        const double lambda = 1.0;
        Bound(Sol, lambda, elem->nc, SplitStep, setZ); 
        /* if(SplitStep < elem->ndim - 1) */
        (*rotate[elem->ndim - 1])(Sol, mpv->Level[0]->rhs, Sbg, buoyS, FORWARD);
    }         
    /* rotate back           
    for(SplitStep = elem->ndim-1; SplitStep > 0; SplitStep--) {
        (*rotate[elem->ndim - 1])(Sol, mpv->Level[0]->rhs, Sbg, buoyS, BACKWARD);
    }
     */
}


/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: boundary.c,v $
 Revision 1.1  1998/03/01 18:43:32  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
