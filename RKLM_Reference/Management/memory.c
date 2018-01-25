/*******************************************************************************
 File:   memory.c
 Author: Rupert
 Date:   ?
 *******************************************************************************/
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "Common.h"
#include "enumerator.h"
#include "kgrid.h"
#include "error.h"
#include "warning.h"
#include "userdata.h"
#include "memory.h"
#include "recovery.h"
#include "variable.h"
#include "Gasdynamics.h"
#include "math_own.h"

/* ================================================================================ */

void update(
			ConsVars* sol, 
			const ConsVars *flux[3], 
            const double* buoyS,
			const VectorField* buoy,
			const ElemSpaceDiscr* elem, 
			const double dt) 
{
	
	extern User_Data ud;
    extern MPV *mpv;
	
    const double grav = ud.gravity_strength[1];
    const double Msq  = ud.Msq;
    const int    impg = ud.implicit_gravity_theta;
    
    double drho, drhoe_total;
  	int ndim = elem->ndim;
    int nsp;
	
  	switch(ndim) {
  		case 1: {
    		int i;
    		const int igx = elem->igx;
    		const int icx = elem->icx;
    		const double lambdax = dt / elem->dx;
    		const ConsVars* f = flux[0];
			
    		for(i = igx; i < icx - igx; i++) {
      			sol->rho[i]  -= lambdax * (f->rho[i + 1]  - f->rho[i]);
      			sol->rhou[i] -= lambdax * (f->rhou[i + 1] - f->rhou[i]) - buoy->x[i];
      			sol->rhov[i] -= lambdax * (f->rhov[i + 1] - f->rhov[i]);
      			sol->rhow[i] -= lambdax * (f->rhow[i + 1] - f->rhow[i]);
      			sol->rhoe[i] -= lambdax * (f->rhoe[i + 1] - f->rhoe[i]);
      			sol->rhoY[i] -= lambdax * (f->rhoY[i + 1] - f->rhoY[i]);
                for (nsp = 0; nsp < ud.nspec; nsp++) {
                    sol->rhoX[nsp][i] -= lambdax * (f->rhoX[nsp][i + 1] - f->rhoX[nsp][i]);
                }
                sol->rhoX[nsp][i] += buoyS[i];
    		}
    		break;
  		}
  		case 2: {
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int ifx = elem->ifx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int ify = elem->ify;
			const double lambdax = dt / elem->dx;
			const double lambday = dt / elem->dy;
			const ConsVars* f = flux[0];
			const ConsVars* g = flux[1];
            
            double drhoY, Yinv, Sxc, Sym, Syp, Nsq, ooopNsqsc;

            int i, j, m, n, ox, oy;

			for(j = igy; j < icy - igy; j++) {
                m = j * icx;
                
                Sxc = mpv->HydroState->S0[j];
                Sym = mpv->HydroState_n->S0[j];
                Syp = mpv->HydroState_n->S0[j+1];
                Nsq = - (grav/Msq) * (Syp-Sym) / (Sxc * elem->dy);
                ooopNsqsc = 1.0 / (1.0 + 0.25 * dt * dt * Nsq);
                
				for(i = igx; i < icx - igx; i++) {
                    n  = m + i;
					ox = j * ifx + i;
					oy = i * ify + j;
					                    
					drho        = - lambdax * (f->rho[ox+1]  - f->rho[ox] ) - lambday * (g->rho[oy+1]  - g->rho[oy] );
                    drhoY       = - lambdax * (f->rhoY[ox + 1] - f->rhoY[ox]) - lambday * (g->rhoY[oy + 1] - g->rhoY[oy]);
                    drhoe_total = - lambdax * (f->rhoe[ox+1] - f->rhoe[ox]) - lambday * (g->rhoe[oy+1] - g->rhoe[oy]);

                    sol->rhoY[n] += drhoY;
					sol->rho[n]  += drho; 
					sol->rhoe[n] += drhoe_total;
                    
                    Yinv = sol->rho[n] / sol->rhoY[n];
					
					sol->rhou[n] += -lambdax * (f->rhou[ox + 1] - f->rhou[ox]) - lambday * (g->rhou[oy + 1] - g->rhou[oy]) + Yinv * buoy->x[n];
					sol->rhov[n] += -lambdax * (f->rhov[ox + 1] - f->rhov[ox]) - lambday * (g->rhov[oy + 1] - g->rhov[oy]) + Yinv * buoy->y[n];
					sol->rhow[n] += -lambdax * (f->rhow[ox + 1] - f->rhow[ox]) - lambday * (g->rhow[oy + 1] - g->rhow[oy]);
                    for (nsp = 0; nsp < ud.nspec; nsp++) {
                        sol->rhoX[nsp][n] -= lambdax * (f->rhoX[nsp][ox + 1] - f->rhoX[nsp][ox]) - lambday * (g->rhoX[nsp][oy + 1] - g->rhoX[nsp][oy]);
                    }
                    
                    sol->rhoX[nsp][n] += buoyS[n];
				}
			}
			            
			break;
		}
		case 3: {
			const int igx = elem->igx;
			const int icx = elem->icx;
			const int ifx = elem->ifx;
			const int igy = elem->igy;
			const int icy = elem->icy;
			const int ify = elem->ify;
			const int igz = elem->igz;
			const int icz = elem->icz;
			const int ifz = elem->ifz;
			const int ifxicy = ifx * icy;
			const int ifyicz = ify * icz;
			const int ifzicx = ifz * icx;
			const double lambdax = dt / elem->dx;
			const double lambday = dt / elem->dy;
			const double lambdaz = dt / elem->dz;
			const ConsVars* f = flux[0];
			const ConsVars* g = flux[1];
			const ConsVars* h = flux[2];
			
            double drhoYS, dYS, Sc, Syp, Sym, rho_old, rhoY_old, rhoYS_old;

            int i, j, k, l, m, n, ox, oy, oz;

            for(k = igz; k < icz - igz; k++) {l = k * icx * icy;
				for(j = igy; j < icy - igy; j++) {m = l + j * icx;
                    
                    Sc  = mpv->HydroState->S0[j];
                    Sym = 0.5*(mpv->HydroState->S0[j]+mpv->HydroState->S0[j-1]);
                    Syp = 0.5*(mpv->HydroState->S0[j]+mpv->HydroState->S0[j+1]);
                    
					for(i = igx; i < icx - igx; i++) {n = m + i;
                        ox = k * ifxicy + j * ifx + i;
                        oy = i * ifyicz + k * ify + j;
                        oz = j * ifzicx + i * ifz + k;

                        /* Buoyancy contribution due to update of theta-fluctuations - preparation */
                        rho_old   = sol->rho[n];
                        rhoY_old  = sol->rhoY[n];
                        rhoYS_old = rho_old - rhoY_old*Sc;
                        drhoYS    = - lambdax * ((f->rho[ox+1] - f->rhoY[ox+1]*Sc) - (f->rho[ox] - f->rhoY[ox]*Sc) ) 
                                    - lambday * ((g->rho[oy+1] - g->rhoY[oy+1]*Syp) - (g->rho[oy] - g->rhoY[oy]*Sym) ) 
                                    - lambdaz * ((h->rho[oz+1] - h->rhoY[oz+1]*Sc) - (h->rho[oz] - h->rhoY[oz]*Sc) );

						/* assuming non-moving grid and stationary geopotential */
						sol->rho[n] -= lambdax * ( f->rho[ox + 1] -  f->rho[ox]);
						sol->rhou[n] -= lambdax * (f->rhou[ox + 1] - f->rhou[ox]) + buoy->x[n];
						sol->rhov[n] -= lambdax * (f->rhov[ox + 1] - f->rhov[ox]);
						sol->rhow[n] -= lambdax * (f->rhow[ox + 1] - f->rhow[ox]);
						sol->rhoe[n] -= lambdax * (f->rhoe[ox + 1] - f->rhoe[ox]);
						sol->rhoY[n] -= lambdax * (f->rhoY[ox + 1] - f->rhoY[ox]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            sol->rhoX[nsp][n] -= lambdax * (f->rhoX[nsp][ox + 1] - f->rhoX[nsp][ox]);
                        }
						/* assuming non-moving grid and stationary geopotential */
						sol->rho[n] -= lambday * ( g->rho[oy + 1] -  g->rho[oy]);
						sol->rhou[n] -= lambday * (g->rhou[oy + 1] - g->rhou[oy]);
						sol->rhov[n] -= lambday * (g->rhov[oy + 1] - g->rhov[oy]) + buoy->y[n];
						sol->rhow[n] -= lambday * (g->rhow[oy + 1] - g->rhow[oy]);
						sol->rhoe[n] -= lambday * (g->rhoe[oy + 1] - g->rhoe[oy]);
						sol->rhoY[n] -= lambday * (g->rhoY[oy + 1] - g->rhoY[oy]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            sol->rhoX[nsp][n] -= lambday * (g->rhoX[nsp][oy + 1] - g->rhoX[nsp][oy]);
                        }

						/* assuming non-moving grid and stationary geopotential */
						sol->rho[n] -= lambdaz * ( h->rho[oz + 1] -  h->rho[oz]);
						sol->rhou[n] -= lambdaz * (h->rhou[oz + 1] - h->rhou[oz]);
						sol->rhov[n] -= lambdaz * (h->rhov[oz + 1] - h->rhov[oz]);
						sol->rhow[n] -= lambdaz * (h->rhow[oz + 1] - h->rhow[oz]) + buoy->z[n];
						sol->rhoe[n] -= lambdaz * (h->rhoe[oz + 1] - h->rhoe[oz]);
						sol->rhoY[n] -= lambdaz * (h->rhoY[oz + 1] - h->rhoY[oz]);
                        for (nsp = 0; nsp < ud.nspec; nsp++) {
                            sol->rhoX[nsp][n] -= lambdaz * (h->rhoX[nsp][oz + 1] - h->rhoX[nsp][oz]);
                        }
                        
                        sol->rhoX[BUOY][n] += buoyS[n];
					}
				} 
			}  
			break;
		}
		default: ERROR("ndim not in {1,2,3}");
	}
}

/* ================================================================================ */

void rotate2D(ConsVars* Sol, double* rhs, double *Sbg, double* buoyS, const enum Direction direction) 
{
	
	/* User data */
	extern User_Data ud;
	
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	extern NodeSpaceDiscr* node;
	
	/* Arrays */
	extern ConsVars* flux[3];             
	extern double *W0;   
	
	double *phelp;
	double delta;
	int i, nsp;
    
	int icx = elem->icx;
	int icy = elem->icy;
	int nc  = elem->nc;
	
	/* rotation of solution */
	flip2D(Sol->rho,  icx, icy, nc, W0);   
	flip2D(Sol->rhou, icx, icy, nc, W0);       
	flip2D(Sol->rhov, icx, icy, nc, W0);       
	flip2D(Sol->rhow, icx, icy, nc, W0);       
	flip2D(Sol->rhoe, icx, icy, nc, W0);       
	flip2D(Sol->rhoY, icx, icy, nc, W0);       
    for (nsp = 0; nsp < ud.nspec; nsp++) {
        flip2D(Sol->rhoX[nsp], icx, icy, nc, W0);
    }
    for (nsp = 0; nsp < ud.naux; nsp++) {
        flip2D(Sol->rhoZ[nsp], icx, icy, nc, W0);
    }
	flip2D(Sol->geopot, icx, icy, nc, W0); 
	
	flip2D(rhs, icx, icy, nc, W0);
    flip2D(Sbg, icx, icy, nc, W0);
    flip2D(buoyS, icx, icy, nc, W0);
	
	/* rotation of grid parameters */
	/* control volumes (elem) */
    i         = elem->nfx;
    elem->nfx = elem->nfy;
    elem->nfy = i;

    i          = elem->igx;
	elem->igx  = elem->igy;
	elem->igy  = i;

    elem->ig[0] = elem->igx;
    elem->ig[1] = elem->igy;

    i          = elem->icx;
    elem->icx  = elem->icy;
    elem->icy  = i;

    elem->ic[0] = elem->icx;
    elem->ic[1] = elem->icy;

    i         = elem->ifx;
	elem->ifx = elem->ify;
	elem->ify = i;

    elem->stride[0] = 1;
    elem->stride[1] = elem->icx;

    delta      = elem->dx;
	elem->dx   = elem->dy;
	elem->dy   = delta;
	
	elem->dxyz[0] = elem->dx;
	elem->dxyz[1] = elem->dy;
	
	phelp   = elem->x;
	elem->x = elem->y;
	elem->y = phelp;
    
    
    /* vertices (node) */
    i         = node->nfx;
    node->nfx = node->nfy;
    node->nfy = i;
    
    i          = node->igx;
    node->igx  = node->igy;
    node->igy  = i;
    
    node->ig[0] = node->igx;
    node->ig[1] = node->igy;
    
    i          = node->icx;
    node->icx  = node->icy;
    node->icy  = i;

    node->ic[0] = node->icx;
    node->ic[1] = node->icy;
	
	i         = node->ifx;
	node->ifx = node->ify;
	node->ify = i;
	
    node->stride[0] = 1;
    node->stride[1] = node->icx;

    delta      = node->dx;
	node->dx   = node->dy;
	node->dy   = delta;

    node->dxyz[0] = node->dx;
    node->dxyz[1] = node->dy;
    
	phelp   = node->x;
	node->x = node->y;
	node->y = phelp;
	    
	/* rotation of momenta and velocities */
	phelp      = Sol->rhou;
	Sol->rhou  = Sol->rhov;
	Sol->rhov  = phelp;
	
	{
		phelp = flux[1]->rhou;
		flux[1]->rhou = flux[1]->rhov;
		flux[1]->rhov = phelp;
	}
}

/* ================================================================================ */

void rotate3D(ConsVars* Sol, double *rhs, double *Sbg, double* buoyS, const enum Direction direction) {
	
	/* User data */
	extern User_Data ud;
	
	/* Grid and space discretization */
	extern ElemSpaceDiscr* elem;
	extern NodeSpaceDiscr* node;
	
	/* Variables */
	extern ConsVars* flux[3];   
	extern double *W0; 
	
	double *phelp;
	double delta;
	int i;
	
	if(direction == FORWARD) {
		
		const int nc  = elem->nc; 
		const int icx = elem->icx;
		const int icy = elem->icy;
		const int icz = elem->icz;
		
		/* rotation of solution arrays */
		flip3D_f( Sol->rho,  icx, icy, icz, nc, W0 );
		flip3D_f( Sol->rhou, icx, icy, icz, nc, W0 );       
		flip3D_f( Sol->rhov, icx, icy, icz, nc, W0 );       
		flip3D_f( Sol->rhow, icx, icy, icz, nc, W0 );       
		flip3D_f( Sol->rhoe, icx, icy, icz, nc, W0 );       
		flip3D_f( Sol->rhoY, icx, icy, icz, nc, W0 );       
        for (int nsp = 0; nsp < ud.nspec; nsp++) {
            flip3D_f(Sol->rhoX[nsp], icx, icy, icz, nc, W0);
        }
        for (int nsp = 0; nsp < ud.naux; nsp++) {
            flip3D_f(Sol->rhoZ[nsp], icx, icy, icz, nc, W0);
        }
		flip3D_f( Sol->geopot, icx, icy, icz, nc, W0 );       
		
		flip3D_f( rhs, icx, icy, icz, nc, W0 );       
        flip3D_f( Sbg, icx, icy, icz, nc, W0 );       
        flip3D_f( buoyS, icx, icy, icz, nc, W0 );       
		
		/* rotation of grid parameters */
		/* control volumes (elem) */
        i         = elem->nfy;
        elem->nfy = elem->nfz;
        elem->nfz = elem->nfx;
        elem->nfx = i;
        
        i         = elem->igy;
        elem->igy = elem->igz;
        elem->igz = elem->igx;
        elem->igx = i;

        elem->ig[0] = elem->igx;
        elem->ig[1] = elem->igy;
        elem->ig[2] = elem->igz;

		i         = elem->icy;
		elem->icy = elem->icz;
		elem->icz = elem->icx;
		elem->icx = i;
		
        elem->ic[0] = elem->icx;
        elem->ic[1] = elem->icy;
        elem->ic[2] = elem->icz;

        i         = elem->ify;
		elem->ify = elem->ifz;
		elem->ifz = elem->ifx;
		elem->ifx = i;

		i               = elem->stride[1];
		elem->stride[1] = elem->stride[2];
		elem->stride[2] = elem->stride[0];
		elem->stride[0] = i;
        
        delta     = elem->dy;
        elem->dy  = elem->dz;
        elem->dz  = elem->dx;
        elem->dx  = delta;
        
        elem->dxyz[0] = elem->dx;
        elem->dxyz[1] = elem->dy;
        elem->dxyz[2] = elem->dz;
        
        phelp   = elem->y;
        elem->y = elem->z;
        elem->z = elem->x;
        elem->x = phelp;
        		
        
		/* vertices (node) */
        i         = node->nfy;
        node->nfy = node->nfz;
        node->nfz = node->nfx;
        node->nfx = i;

        i         = node->igy;
        node->igy = node->igz;
        node->igz = node->igx;
        node->igx = i;

        node->ig[0] = node->igx;
        node->ig[1] = node->igy;
        node->ig[2] = node->igz;

        i         = node->icy;
		node->icy = node->icz;
		node->icz = node->icx;
		node->icx = i;
		
        node->ic[0] = node->icx;
        node->ic[1] = node->icy;
        node->ic[2] = node->icz;
        
        i         = node->ify;
		node->ify = node->ifz;
		node->ifz = node->ifx;
		node->ifx = i;
        
        i               = node->stride[1];
        node->stride[1] = node->stride[2];
        node->stride[2] = node->stride[0];
        node->stride[0] = i;        

        delta     = node->dy;
		node->dy  = node->dz;
		node->dz  = node->dx;
		node->dx  = delta;

        node->dxyz[0] = node->dx;
        node->dxyz[1] = node->dy;
        node->dxyz[2] = node->dz;

		phelp   = node->y;
		node->y = node->z;
		node->z = node->x;
		node->x = phelp;
		
		/* rotation of momenta and velocities */
		phelp     = Sol->rhov;
		Sol->rhov = Sol->rhow;
		Sol->rhow = Sol->rhou;
		Sol->rhou = phelp;
		
		{
			phelp = flux[0]->rhov;
			flux[0]->rhov = flux[0]->rhow;
			flux[0]->rhow = flux[0]->rhou;
			flux[0]->rhou = phelp;
			phelp = flux[1]->rhov;
			flux[1]->rhov = flux[1]->rhow;
			flux[1]->rhow = flux[1]->rhou;
			flux[1]->rhou = phelp;
			phelp = flux[2]->rhov;
			flux[2]->rhov = flux[2]->rhow;
			flux[2]->rhow = flux[2]->rhou;
			flux[2]->rhou = phelp;
		}
	}
	else if(direction == BACKWARD) {
		
		int nc, icx, icy, icz;
		
		/* rotation of grid increments */
		/* control volumes (elem) */
        i         = elem->nfz;
        elem->nfz = elem->nfy;
        elem->nfy = elem->nfx;
        elem->nfx = i;
        
        i         = elem->igz;
        elem->igz = elem->igy;
        elem->igy = elem->igx;
        elem->igx = i;
        
        elem->ig[0] = elem->igx;
        elem->ig[1] = elem->igy;
        elem->ig[2] = elem->igz;

        i         = elem->icz;
		elem->icz = elem->icy;
		elem->icy = elem->icx;
		elem->icx = i;
        
        elem->ic[0] = elem->icx;
        elem->ic[1] = elem->icy;
        elem->ic[2] = elem->icz;
		
		i         = elem->ifz;
		elem->ifz = elem->ify;
		elem->ify = elem->ifx;
		elem->ifx = i;

        i               = elem->stride[2];
        elem->stride[2] = elem->stride[1];
        elem->stride[1] = elem->stride[0];
        elem->stride[0] = i;        

        delta    = elem->dz;
		elem->dz = elem->dy;
		elem->dy = elem->dx; 
		elem->dx = delta;
		
		elem->dxyz[0] = elem->dx;
		elem->dxyz[1] = elem->dy;
		elem->dxyz[2] = elem->dz;
		
		phelp   = elem->z;
		elem->z = elem->y;
		elem->y = elem->x;
		elem->x = phelp;
        

        /* vertices (node) */
        i         = node->nfz;
        node->nfz = node->nfy;
        node->nfy = node->nfx;
        node->nfx = i;
        
        i         = node->igz;
        node->igz = node->igy;
        node->igy = node->igx;
        node->igx = i;
        
        node->ig[0] = node->igx;
        node->ig[1] = node->igy;
        node->ig[2] = node->igz;

        i         = node->icz;
		node->icz = node->icy;
		node->icy = node->icx;
		node->icx = i;
		
        node->ic[0] = node->icx;
        node->ic[1] = node->icy;
        node->ic[2] = node->icz;
        
		i         = node->ifz;
		node->ifz = node->ify;
		node->ify = node->ifx;
		node->ifx = i;

        i               = node->stride[2];
        node->stride[2] = node->stride[1];
        node->stride[1] = node->stride[0];
        node->stride[0] = i;        

        delta    = node->dz;
		node->dz = node->dy;
		node->dy = node->dx; 
		node->dx = delta;
		
        node->dxyz[0] = node->dx;
        node->dxyz[1] = node->dy;
        node->dxyz[2] = node->dz;

        phelp   = node->z;
		node->z = node->y;
		node->y = node->x;
		node->x = phelp;
                
		
		/* rotation of momenta and velocities */
		phelp     = Sol->rhow;
		Sol->rhow = Sol->rhov;
		Sol->rhov = Sol->rhou;
		Sol->rhou = phelp;
		
		{
			phelp = flux[0]->rhow;
			flux[0]->rhow = flux[0]->rhov;
			flux[0]->rhov = flux[0]->rhou;
			flux[0]->rhou = phelp;
			phelp = flux[1]->rhow;
			flux[1]->rhow = flux[1]->rhov;
			flux[1]->rhov = flux[1]->rhou;
			flux[1]->rhou = phelp;
			phelp = flux[2]->rhow;
			flux[2]->rhow = flux[2]->rhov;
			flux[2]->rhov = flux[2]->rhou;
			flux[2]->rhou = phelp;
		}
		
		nc  = elem->nc; 
		icx = elem->icx;
		icy = elem->icy;
		icz = elem->icz;
		
		/* rotation of solution arrays */
		flip3D_b( Sol->rho,  icx, icy, icz, nc, W0 );       
		flip3D_b( Sol->rhou, icx, icy, icz, nc, W0 );       
		flip3D_b( Sol->rhov, icx, icy, icz, nc, W0 );       
		flip3D_b( Sol->rhow, icx, icy, icz, nc, W0 );       
		flip3D_b( Sol->rhoe, icx, icy, icz, nc, W0 );       
		flip3D_b( Sol->rhoY, icx, icy, icz, nc, W0 );       
        for (int nsp = 0; nsp < ud.nspec; nsp++) {
            flip3D_b(Sol->rhoX[nsp], icx, icy, icz, nc, W0);
        }
        for (int nsp = 0; nsp < ud.naux; nsp++) {
            flip3D_b(Sol->rhoZ[nsp], icx, icy, icz, nc, W0);
        }
		flip3D_b( Sol->geopot, icx, icy, icz, nc, W0 );       
		
		flip3D_b( rhs, icx, icy, icz, nc, W0 );       
        flip3D_b( Sbg, icx, icy, icz, nc, W0 );       
        flip3D_b( buoyS, icx, icy, icz, nc, W0 );       
		
	}
	else {
		ERROR("which stupid direction shall i rotate the stuff, eehhh ?");
	}
}  

/* ================================================================================ */

void flip3D_f( double *f, int ix, int iy, int iz, int n, double *W ) {
	
	double *pf, *paux; 
	int i, j, k;
	
	for(i = 0, pf = f, paux = W; i < n; i++, pf++, paux++)
		*paux = *pf;
	
	for(k = 0; k < iz; k++)
		for(j = 0; j < iy; j++)
			for(i = 0; i < ix; i++)
				f[i * iy * iz + k * iy + j] = W[k * ix * iy + j * ix + i];
}

/* ================================================================================ */

void flip3D_b( double *f, int ix, int iy, int iz, int n, double *W ) {
	
	double *pf, *paux; 
	int i, j, k;
	
	for(i = 0, pf = f, paux = W; i < n; i++, pf++, paux++)
		*paux = *pf;
	
	for(k = 0; k < iz; k++)
		for(j = 0; j < iy; j++)
			for(i = 0; i < ix; i++)
				f[k * ix * iy + j * ix + i] = W[i * iy * iz + k * iy + j];
}


/* ================================================================================ */

void flip2D(double *f, int ix, int iy, int n, double *W ) {
	
	int i, j;
	
	memcpy(W, f, n * sizeof(double));
	
	for(j = 0; j < iy; j++)
		for(i = 0; i < ix; i++)
			f[i * iy + j] = W[j * ix + i];
}

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: memory.c,v $
 Revision 1.2  1998/03/07 09:56:46  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:35  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
