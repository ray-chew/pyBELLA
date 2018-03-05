/*******************************************************************************
 File:   memory.c
 Author: Rupert
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

void rotate2D(ConsVars* Sol, const enum Direction direction) 
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

void rotate3D(ConsVars* Sol, const enum Direction direction) {
	
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
	}
	else {
		ERROR("to which direction shall i rotate the fields?");
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

/* LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL */
