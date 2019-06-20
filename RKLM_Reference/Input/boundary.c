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
#include "Eos.h"
#include "thermodynamic.h"
#include "math_own.h"
#include "memory.h"

#if OUTPUT_SUBSTEPS 
#include "io.h"
#endif


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

double slanted_wall_slope(double x);
double slanted_wall_slope_3D(double x, double y);


/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/

void initialize_bdry(
					 const ElemSpaceDiscr* elem)
{
	const int icx = elem->icx;
    const int icz = elem->icz;
	const int igx = elem->igx;
    const int igz = elem->igz;

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
            bdry->wall_rhoYflux = (double*)malloc(icx*sizeof(double));
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
            /* note that the y-direction is my "vertical" */
            bdry->wall_rhoYflux = (double*)malloc(icx*icz*sizeof(double));
            bdry->wall_slope    = (double*)malloc(icx*icz*sizeof(double));
            bdry->wall_relative_slope = (double*)malloc(icx*icz*sizeof(double));
            
            slope_sum = 0.0;
            for (int k = igz; k < icz-igz; k++) {
                int nk = k*icx;
                for (int i = igx; i < icx-igx; i++) {
                    int nik = nk + i;
                    double slope = slanted_wall_slope_3D(elem->x[i], elem->z[k]);
                    bdry->wall_slope[nik] = slope;
                    slope_sum += fabs(slope);
                }
            }
            
            if(slope_sum <= sqrt(DBL_EPSILON)) {
                for (int k = igz; k < icz-igz; k++) {
                    int nk = k*icx;
                    for (int i = igx; i < icx-igx; i++) {
                        int nik = nk + i;
                        bdry->wall_relative_slope[nik] = 0.0;
                    }
                }
            }
            else {
                slope_sum_inv = 1.0/slope_sum;
                for (int k = igz; k < icz-igz; k++) {
                    int nk = k*icx;
                    for (int i = igx; i < icx-igx; i++) {
                        int nik = nk + i;
                        bdry->wall_relative_slope[nik] = slope_sum_inv * fabs(bdry->wall_slope[nik]);
                    }
                }
                
                /*	*/
                slope_sum = 0.0;
                for (int k = igz; k < icz-igz; k++) {
                    int nk = k*icx;
                    for (int i = igx; i < icx-igx; i++) {
                        int nik = nk + i;
                        slope_sum += bdry->wall_relative_slope[nik];
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
    free(bdry->wall_rhoYflux);
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
		   const int SplitStep) 
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
					double rhoYu_wall, rhoYu_image;
					
                    /* copy wall-tangential velocities and scalars */
                    v = Sol->rhov[nsource] / Sol->rho[nsource]; 
                    w = Sol->rhow[nsource] / Sol->rho[nsource];
                    for (nsp = 0; nsp < NSPEC; nsp++) {
                        X[nsp] = Sol->rhoX[nsp][nsource] / Sol->rho[nsource];
                    }
                    
                    /* mirror wall-normal velocity relative to prescribed boundary mass flux */
					rhoYu_wall  = bdry->wall_rhoYflux[j];
                    rhoYu_image = 2.0*rhoYu_wall - Sol->rhou[nsource]*Sol->rhoY[nsource]/Sol->rho[nsource];
					
					{
						/* double S = (2.0*Sol->rho[nimage+1]/Sol->rhoY[nimage+1] - Sol->rho[nimage+2]/Sol->rhoY[nimage+2]); */
                        double S = 1.0;
                        
                        if (ud.bottom_theta_bc == ZERO_ORDER_EXTRAPOL) {
                            /* for the Straka test */
                            S    = 1.0/Y_last; 
                        } else {
                            /* for all other tests so far */
                            int iimage  = i;
                            S    = 1./stratification(elem->x[iimage]); 
                        }
                        double dpi  = (th.Gamma*g) * 0.5*dh*(1.0/Y_last + S);
                        double rhoY = (compressible == 1 ? pow(pow(Sol->rhoY[nlast],th.gm1) + dpi, th.gm1inv) : mpv->HydroState->rhoY0[i]);
                        double rho  = rhoY * S;
                        double p    = pow(rhoY, th.gamm);
                        double u    = rhoYu_image/rhoY;

                        Sol->rho[nimage]  = rho;
						Sol->rhou[nimage] = rho*u;
						Sol->rhov[nimage] = rho*v;
						Sol->rhow[nimage] = rho*w;
						Sol->rhoe[nimage] = rhoe(rho, u, v, w, p);
						Sol->rhoY[nimage] = rhoY;				  /* should probably be adjusted not to take HydroState values*/
                        for (nsp = 0; nsp < NSPEC; nsp++) {
                            Sol->rhoX[nsp][nimage] = rho * X[nsp];
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
						Sol->rhoe[nimage] = rhoe(rho, u, v, w, p);
						Sol->rhoY[nimage] = rhoY;       /* should probably be adjusted not to take HydroState values*/
                        for (nsp = 0; nsp < NSPEC; nsp++) {
                            Sol->rhoX[nsp][nimage] = rho * X[nsp];
                        }
					}
                }
            }
        }
    }
}

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/
void set_wall_rhoYflux(
					   BDRY* bdry, 
					   const ConsVars* Sol0, 
                       const MPV* mpv,
					   const ElemSpaceDiscr* elem)
{
    extern User_Data ud;
    
	const int icx = elem->icx;
	const int icy = elem->icy;
    const int icz = elem->icz;
	const int igx = elem->igx;
    const int igy = elem->igy;
    const int igz = elem->igz;

    double wall_mf, wall_flux_balance;
	int i;
	
    switch (elem->ndim) {
        case 1: {
            printf("\nwall flux in 1D makes no sense\n");
            exit(15);
            break;                
        }
            
        case 2: {
            int is_x_periodic = 0;            
            if(ud.bdrytype_min[0] == PERIODIC) is_x_periodic = 1;

            /* first guess for wall mass fluxes */
            const int nstart = elem->igy*elem->icx;

            wall_flux_balance = 0.0;
            for (i = igx; i < icx-igx; i++) {
                int n = nstart + i;
                wall_mf = wall_rhoYflux(elem->x[i], elem->y[i], Sol0->rhou[n]/Sol0->rho[n], Sol0->rhov[n]/Sol0->rho[n], mpv->HydroState_n->rhoY0[igy]);
                bdry->wall_rhoYflux[i] = wall_mf;
                wall_flux_balance += wall_mf;
            }
            
            /* correction for zero net flux */
            for (i = igx; i < icx-igx; i++) {
                bdry->wall_rhoYflux[i] -= wall_flux_balance * bdry->wall_relative_slope[i];
            }

            for(i=0; i<elem->igx; i++) {
                bdry->wall_rhoYflux[i] = bdry->wall_rhoYflux[icx-1-i] = 0.0; 
            }
            if (is_x_periodic) {
                for(i=0; i<elem->igx; i++) {
                    bdry->wall_rhoYflux[i] = bdry->wall_rhoYflux[icx-2*igx+i]; 
                    bdry->wall_rhoYflux[icx-igx+i] = bdry->wall_rhoYflux[igx+i]; 
                }
            }
            
            /*	*/
            {
                double flux_sum = 0.0;
                for (i = elem->igx; i < elem->icx-elem->igx; i++) {
                    flux_sum += bdry->wall_rhoYflux[i];
                }
                
                /* printf("wall flux sum = %e\n", flux_sum); */
            }
            break;
        }
            
        case 3: {
            int is_x_periodic = 0;            
            int is_z_periodic = 0;            
            if(ud.bdrytype_min[0] == PERIODIC) is_x_periodic = 1;
            if(ud.bdrytype_min[2] == PERIODIC) is_z_periodic = 1;

            /* y-direction is vertical; wall flux bdry is the bottom  x-z-surface */
            const int nstart = elem->igy*elem->icx;

            /* first guess for wall mass fluxes before total flux correction (for elliptic solvability cond.) */
            wall_flux_balance = 0.0;
            for (int k = igz; k < icz-igz; k++) {
                int njk = nstart + k*icx*icy;
                int nk  = k*icx;
                for (int i = igx; i < icx-igx; i++) {
                    int nijk = njk + i;
                    int nik  = nk  + i;
                    wall_mf = wall_rhoYflux(elem->x[i], elem->z[k], Sol0->rhou[nijk]/Sol0->rho[nijk], Sol0->rhow[nijk]/Sol0->rho[nijk], mpv->HydroState_n->rhoY0[igy]);
                    bdry->wall_rhoYflux[nik] = wall_mf;
                    wall_flux_balance += wall_mf;
                }
            }
            
            /* set wall flux in dummy cells to zero */
            for(int k=0; k<icz; k++) {
                int nk = k * icx;
                for(int i=0; i<igx; i++) {
                    int nik_left  = nk + i;
                    int nik_right = nk + icx-1-i;
                    bdry->wall_rhoYflux[nik_left] = bdry->wall_rhoYflux[nik_right] = 0.0; 
                }
            }            

            for(int i=0; i<icx; i++) {
                int ni = i;
                for(int k=0; k<igz; k++) {
                    int nik_left  = ni + k*icx;
                    int nik_right = ni + (icz-1-k)*icx;
                    bdry->wall_rhoYflux[nik_left] = bdry->wall_rhoYflux[nik_right] = 0.0; 
                }
            }            
            
            /* correction for zero net flux (needed for elliptic solvability) */
            for (int k = igz; k < icz-igz; k++) {
                int nk = k * icx;
                for (int i = igx; i < icx-igx; i++) {
                    int nik  = nk + i;
                    bdry->wall_rhoYflux[nik] -= wall_flux_balance * bdry->wall_relative_slope[nik];
                }
            }
            
            /* the remaining lines in this routine have not yet been tested */
            for (int k=0; k<icz; k++) {
                int nk = k*icx;
                for(int i=0; i<elem->igx; i++) {
                    int niklt = nk + i;
                    int nikrt = nk + icx-igx+i;
                    bdry->wall_rhoYflux[niklt] = 0.0; 
                    bdry->wall_rhoYflux[nikrt] = 0.0; 
                }
            }
            if (is_x_periodic) {
                for (int k=0; k<icz; k++) {
                    int nk = k*icx;
                    for(int i=0; i<elem->igx; i++) {
                        int niklt = nk + i;
                        int nikrs = nk + icx-2*igx+i;
                        int nikrt = nk + icx-igx+i;
                        int nikls = nk + igx+i;
                        bdry->wall_rhoYflux[niklt] = bdry->wall_rhoYflux[nikrs]; 
                        bdry->wall_rhoYflux[nikrt] = bdry->wall_rhoYflux[nikls]; 
                    }
                }
            } 
            
            for (int i=0; i<icx; i++) {
                int ni = i;
                for(int k=0; k<elem->igz; k++) {
                    int niklt = ni + k*icx;
                    int nikrt = ni + (icz-igz+k)*icx;
                    bdry->wall_rhoYflux[niklt] = 0.0; 
                    bdry->wall_rhoYflux[nikrt] = 0.0; 
                }
            }
            if (is_z_periodic) {
                for (int i=0; i<icx; i++) {
                    int ni = i;
                    for(int k=0; k<elem->igz; k++) {
                        int niklt = ni + k*icx;
                        int nikrs = ni + (icz-2*igz+k)*icx;
                        int nikrt = ni + (icz-igz+k)*icx;
                        int nikls = ni + (igz+k)*icx;
                        bdry->wall_rhoYflux[niklt] = bdry->wall_rhoYflux[nikrs]; 
                        bdry->wall_rhoYflux[nikrt] = bdry->wall_rhoYflux[nikls]; 
                    }
                }
            } 

            /*	*/
            {
                double flux_sum = 0.0;
                for (int k = igz; k < icz-igz; k++) {
                    int nk = k * icx;
                    for (int i = igx; i < icx-igx; i++) {
                        int nik  = nk + i;
                        flux_sum += bdry->wall_rhoYflux[nik];
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
	double u_outer   = velo_background(t);
	
	nijk    = njk + elem->igx-1 - i;
	
	Sol->rho[nijk]    = rho_outer;
	Sol->rhou[nijk]   = rho_outer * u_outer;
	Sol->rhov[nijk]   = 0.0;
	Sol->rhow[nijk]   = 0.0;
	Sol->rhoe[nijk]   = rhoe(rho_outer, u_outer, 0.0, 0.0, p_outer);
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
	double u_outer   = velo_background(t);
	
	nijk = njk + i;
	
	Sol->rho[nijk]    = rho_outer;
	Sol->rhou[nijk]   = rho_outer * u_outer;
	Sol->rhov[nijk]   = 0.0;
	Sol->rhow[nijk]   = 0.0;
	Sol->rhoe[nijk]   = rhoe(rho_outer, u_outer, 0.0, 0.0, p_outer);
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
    extern double t;
    
	double hill_height  = ud.hill_height;
	double length_scale_inv = 1.0/ud.hill_length_scale; 
	double x_sc = x*length_scale_inv;
	
    if (ud.hill_shape == SCHLUTOW) {
        
        /* Topography for Mark Schlutow's stationary WKB waves */
        double kx = 2.0*2.0*PI/(ud.xmax-ud.xmin);
        double kz = sqrt(ud.Nsq/ud.wind_speed/ud.wind_speed + kx*kx);
        double q  = 0.25;

        /* scaled height and slope */
        double xi = kx*x;
        double y  = 0.0;
        double yp;
        
        for (int i=0; i<10; i++) {
            y = q*cos(xi+y);
        }
        yp = - q*sin(xi+y)/(1+q*sin(xi+y));
        
        /* unscaled slope */
        yp = kx*yp/kz;
        return(yp);
    } else {
        return(- hill_height * 2.0*x_sc / ((1.0 + x_sc*x_sc)*(1.0 + x_sc*x_sc)) * length_scale_inv);        
    }
	
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

/* ========================================================================= */

double wall_rhoYflux(const double x, 
                     const double z, 
                     const double wind_speed_x, 
                     const double wind_speed_z, 
                     const double rhoY0){
	extern User_Data ud;
    	
    if (ud.hill_shape == SCHLUTOW) {
        return (slanted_wall_slope(x) * wind_speed_x * rhoY0);
    } else {
#if 1
        return (slanted_wall_slope(x) * wind_speed_x * rhoY0);
#else
        double hill_height  = ud.hill_height;
        double length_scale_inv = 1.0/ud.hill_length_scale; 
        double x_sc = x*length_scale_inv;
        
        return((- wind_speed_x * hill_height * 2.0*x_sc / ((1.0 + x_sc*x_sc)*(1.0 + x_sc*x_sc)) * length_scale_inv
                - wind_speed_z * 0.0) * rhoY0);
#endif
    }
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
		
    int i, nfull, ifull, jfull, nsp; 
    
	switch (SplitStep) {
		case 1: {
			if(ud.bdrytype_min[SplitStep] == WALL) {
				for(i=elem->igx; i<nmax-elem->igx+1; i++){
					nfull = kcache*njump + i;
					jfull = nfull/elem->icx;
					ifull = nfull%elem->icx;
					if(ifull==elem->igx) {
                        double rhou_wall = bdry->wall_rhoYflux[jfull];
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
						double S2  = mpv->HydroState->Y0[jfull];
						
						Rights->rho[i]  = Lefts->rho[i-1]  = rho;
						Rights->rhou[i] = Lefts->rhou[i-1] = rho * ud.wind_speed;
						Rights->rhov[i] = Lefts->rhov[i-1] = 0.0;
						Rights->rhow[i] = Lefts->rhow[i-1] = 0.0;
						Rights->rhoe[i] = Lefts->rhoe[i-1] = rhoe(rho, ud.wind_speed, 0.0, 0.0, p);
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
						double S2  = mpv->HydroState->Y0[jfull];
						
						Lefts->rho[i-1]  = Rights->rho[i]  = rho;
						Lefts->rhou[i-1] = Rights->rhou[i] = rho * ud.wind_speed;
						Lefts->rhov[i-1] = Rights->rhov[i] = 0.0;
						Lefts->rhow[i-1] = Rights->rhow[i] = 0.0;
						Lefts->rhoe[i-1] = Rights->rhoe[i] = rhoe(rho, ud.wind_speed, 0.0, 0.0, p);
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

static void (*rotate[])(ConsVars* Sol, const enum Direction dir) = {NULL, rotate2D, rotate3D};

void Set_Explicit_Boundary_Data(
                                ConsVars* Sol,
                                const ElemSpaceDiscr* elem) 
{
    int SplitStep;
    
    for(SplitStep = 0; SplitStep < elem->ndim; SplitStep++) { 
        const double lambda = 1.0;
        Bound(Sol, lambda, elem->nc, SplitStep); 
        (*rotate[elem->ndim - 1])(Sol, FORWARD);
    }   
    
#if OUTPUT_SUBSTEPS /* 5 */
    extern User_Data ud;
    extern NodeSpaceDiscr* node;
    extern int step;
    if (step >= OUTPUT_SUBSTEPS - 1) {
        putout(Sol, ud.file_name, "Sol", elem, node, 1);
    }
#endif

}

/* ============================================================================= */

void set_ghostcells_p2(
                       double* p,
                       const ElemSpaceDiscr* elem, 
                       const int ig) {
    
    extern User_Data ud;
    
    /* I am implementing a simple, rudimentary version just for homogeneous
     Neumann and periodic conditions
     */

    const int icx = elem->icx;
    const int igx = elem->igx;
    const int icy = elem->icy;
    const int igy = elem->igy;
    const int icz = elem->icz;
    const int igz = elem->igz;
    
    const int xperiodic = (ud.bdrytype_min[0] == PERIODIC ? 1 : 0);
    const int yperiodic = (ud.bdrytype_min[1] == PERIODIC ? 1 : 0);
    const int zperiodic = (ud.bdrytype_min[2] == PERIODIC ? 1 : 0);
    
    /* x-direction */
    for (int k=0; k<icz; k++) {
        int lc = k*icx*icy;
        for (int j=0; j<icy; j++) {
            int mc  = lc + j*icx;
            for (int i=igx-1; i>=igx-ig; i--) {
                int nc0_obj = mc + i;
                int nc0_src = mc + xperiodic*(icx-4+i) + (1-xperiodic)*(igx+1-i);
                int nc1_obj = mc + icx-1-i;
                int nc1_src = mc + xperiodic*(igx+1-i) + (1-xperiodic)*(icx-4+i);
                
                p[nc0_obj] = p[nc0_src];
                p[nc1_obj] = p[nc1_src];
            }
        }
    } 
    
    /* y-direction */
    if (elem->ndim > 1) {
        for (int i=0; i<icx; i++) {
            int lc = i;
            for (int k=0; k<icz; k++) {
                int mc  = lc + k*icx*icy;
                for (int j=igy-1; j>=igy-ig; j--) {
                    int nc0_obj = mc + j*icx;
                    int nc0_src = mc + (yperiodic*(icy-4+j) + (1-yperiodic)*(igy+1-j))*icx;
                    int nc1_obj = mc + (icy-1-j)*icx;
                    int nc1_src = mc + (yperiodic*(igy+1-j) + (1-yperiodic)*(icy-4+j))*icx;
                    
                    p[nc0_obj] = p[nc0_src];
                    p[nc1_obj] = p[nc1_src];
                }
            }
        } 
    }

    /* z-direction */
    if (elem->ndim > 2) {
        for (int j=0; j<icy; j++) {
            int lc = j*icx;
            for (int i=0; i<icx; i++) {
                int mc  = lc + i;
                for (int k=igz-1; k>=igz-ig; k--) {
                    int nc0_obj = mc + k*icx*icy;
                    int nc0_src = mc + (zperiodic*(icz-4+k) + (1-zperiodic)*(igz+1-k))*icx*icy;
                    int nc1_obj = mc + (icz-1-k)*icx*icy;
                    int nc1_src = mc + (zperiodic*(igz+1-k) + (1-zperiodic)*(icz-4+k))*icx*icy;
                        
                    p[nc0_obj] = p[nc0_src];
                    p[nc1_obj] = p[nc1_src];
                }
            }
        } 
    }
}

/* ============================================================================= */

void set_ghostnodes_p2(
                       double* p,
                       const NodeSpaceDiscr* node, 
                       const int ig) {
    
    extern User_Data ud;
    
    /* I am implementing a simple, rudimentary version just for homogeneous
     Neumann and periodic conditions
     */
    
    const int icx = node->icx;
    const int igx = node->igx;
    const int icy = node->icy;
    const int igy = node->igy;
    const int icz = node->icz;
    const int igz = node->igz;
    
    const int xperiodic = (ud.bdrytype_min[0] == PERIODIC ? 1 : 0);
    const int yperiodic = (ud.bdrytype_min[1] == PERIODIC ? 1 : 0);
    const int zperiodic = (ud.bdrytype_min[2] == PERIODIC ? 1 : 0);
    
    /* x-direction */
    for (int k=0; k<icz; k++) {
        int lc = k*icx*icy;
        for (int j=0; j<icy; j++) {
            int mc  = lc + j*icx;
            for (int i=igx; i>=igx-ig; i--) {
                int nc0_obj = mc + i;
                int nc0_src = mc + xperiodic*(icx-5+i) + (1-xperiodic)*(igx+2-i);
                int nc1_obj = mc + icx-1-i;
                int nc1_src = mc + xperiodic*(igx+2-i) + (1-xperiodic)*(icx-5+i);
                
                /* for the periodic case it the sequence "first 1 then 0" should count */
                p[nc1_obj] = p[nc1_src];
                p[nc0_obj] = p[nc0_src];
            }
        }
    } 
    
    /* y-direction */
    if (node->ndim > 1) {
        for (int i=0; i<icx; i++) {
            int lc = i;
            for (int k=0; k<icz; k++) {
                int mc  = lc + k*icx*icy;
                for (int j=igy; j>=igy-ig; j--) {
                    int nc0_obj = mc + j*icx;
                    int nc0_src = mc + (yperiodic*(icy-5+j) + (1-yperiodic)*(igy+2-j))*icx;
                    int nc1_obj = mc + (icy-1-j)*icx;
                    int nc1_src = mc + (yperiodic*(igy+2-j) + (1-yperiodic)*(icy-5+j))*icx;
                    
                    p[nc1_obj] = p[nc1_src];
                    p[nc0_obj] = p[nc0_src];
                }
            }
        } 
    }
    
    /* z-direction */
    if (node->ndim > 2) {
        for (int j=0; j<icy; j++) {
            int lc = j*icx;
            for (int i=0; i<icx; i++) {
                int mc  = lc + i;
                for (int k=igz; k>=igz-ig; k--) {
                    int nc0_obj = mc + k*icx*icy;
                    int nc0_src = mc + (zperiodic*(icz-5+k) + (1-zperiodic)*(igz+2-k))*icx*icy;
                    int nc1_obj = mc + (icz-1-k)*icx*icy;
                    int nc1_src = mc + (zperiodic*(igz+2-k) + (1-zperiodic)*(icz-5+k))*icx*icy;
                    
                    p[nc1_obj] = p[nc1_src];
                    p[nc0_obj] = p[nc0_src];
                }
            }
        } 
    }
}
/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: boundary.c,v $
 Revision 1.1  1998/03/01 18:43:32  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
