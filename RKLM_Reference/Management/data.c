/*******************************************************************************
 File:   data.c
 Author: Nicola                        Rupert 
 Date:   Fri Feb 27 10:42:17 CET 1998  Feb. 2004
 *******************************************************************************/
#include <stdlib.h>
#include <string.h>
#include "Common.h"
#include "enumerator.h"
#include "math_own.h"
#include "variable.h"
#include "userdata.h"
#include "data.h"
#include "kgrid.h"
#include "thermodynamic.h"
#include "Gasdynamics.h"
#include "recovery.h"
#include "boundary.h"
#include "SimpleUtilities.h"
#include "mpv.h"

/*------------------------------------------------------------------------------
 Non-constant data
 ------------------------------------------------------------------------------*/

/* User data */
User_Data ud;

/* Time discretization */
int step;
double t;
double dt;

/* Grid and space discretization */
ElemSpaceDiscr* elem;
NodeSpaceDiscr* node;

/* Thermodynamics */
Thermodynamic th;

/* Arrays */
ConsVars* Sol;            /* full size */
ConsVars* Sol0;           /* full size (M < 1.0) */
ConsVars* dSol;           /* full size */ /* TODO: Can I work without full-size dSol arrays? 
                                           -> No: currently still needed in Explicit_step_update */

double* force[3];
ConsVars* flux[3];        /* full size (M < 1.0) */

#ifdef SYMMETRIC_ADVECTION
ConsVars* Sol1;           /* full size (M < 1.0; needed to test x-y-symmetric advection for the travelling vortex) */
#endif

/* Infrastructure for semi-implicit scheme */
MPV* mpv;

States* Solk;             /* cache size */
double* W0;               /* full nG-length double array */
enum Boolean W0_in_use = WRONG;


/*------------------------------------------------------------------------------
 acquire memory
 ------------------------------------------------------------------------------*/

void Data_init() {
	
    extern BDRY* bdry;

	int inx, iny, inz, n_aux;
	double x0, x1, y0, y1, z0, z1;
	enum BdryType left, right, bottom, top, back, front;
	Grid* grid;
	
	/* User data / problem parameters */
	User_Data_init(&ud);
	
    /* allocate memory for code infrastructure */
	inx = ud.inx;
	iny = ud.iny;
	inz = ud.inz;
	x0 = ud.xmin;
	x1 = ud.xmax;
	y0 = ud.ymin;
	y1 = ud.ymax;
	z0 = ud.zmin;
	z1 = ud.zmax;
	left = ud.bdrytype_min[0];
	right = ud.bdrytype_max[0];
	bottom = ud.bdrytype_min[1];
	top = ud.bdrytype_max[1];
	back = ud.bdrytype_min[2];
	front = ud.bdrytype_max[2];
	
	/* Time discretization */
	step = 0;
	t = 0.0;
	
	/* Grid and space discretization */
	grid = Grid_new(inx,iny,inz,x0,x1,y0,y1,z0,z1,left,right,bottom,top,back,front); 
	elem = ElemSpaceDiscr_new(grid);
	node = NodeSpaceDiscr_new(grid);
	Grid_free(grid); 
	
	/* Thermodynamics */
	Thermodynamic_init(&th, &ud);
	
	/* Arrays */
    Sol0    = ConsVars_new(elem->nc);
	Sol  = ConsVars_new(elem->nc);
	dSol = ConsVars_new(elem->nc);
	Solk = States_small_new(3 * ud.ncache / 2);
    
    /* auxiliary space */
    n_aux = node->ifx;
    n_aux *= (node->ndim > 1 ? node->ify : 1);
    n_aux *= (node->ndim > 2 ? node->ifz : 1);
    
    W0  = (double*)malloc((unsigned)(n_aux * sizeof(double)));
    
    for (int idim = 0; idim < elem->ndim; idim++) {
        force[idim]    = (double*)malloc(node->nc*(sizeof(double)));
        for (int nc = 0; nc < node->nc; nc++) force[idim][nc] = 0.0;
    }
    
    flux[0] = ConsVars_new(elem->nfx);
    ConsVars_setzero(flux[0], elem->nfx);
    if(elem->ndim > 1) {
        flux[1] = ConsVars_new(elem->nfy);
        ConsVars_setzero(flux[1], elem->nfy);
    }
    if(elem->ndim > 2) {
        flux[2] = ConsVars_new(elem->nfz);
        ConsVars_setzero(flux[2], elem->nfz);
    }
    
    initialize_bdry(elem);
    
    mpv = MPV_new(elem, node);
            
    Explicit_malloc(3 * ud.ncache / 2);
	recovery_malloc(3 * ud.ncache / 2);  

#ifdef SYMMETRIC_ADVECTION
    Sol1    = ConsVars_new(elem->nc);
    ConsVars_setzero(Sol1, elem->nc);
#endif

    /* set initial data */
    Sol_initial(Sol, Sol0, mpv, bdry, elem, node);
    ConsVars_setzero(dSol, elem->nc);

}


/*------------------------------------------------------------------------------
 release memory
 ------------------------------------------------------------------------------*/

void Data_free() {
	    
    recovery_free();    
    Explicit_free();

    MPV_free(mpv, elem);

    close_bdry();

    for (int idim = 0; idim < elem->ndim; idim++) {
        ConsVars_free(flux[idim]);
        free(force[idim]);
    }

    free(W0);

    States_small_free(Solk);
    ConsVars_free(dSol);
    ConsVars_free(Sol);
#ifdef SYMMETRIC_ADVECTION
    ConsVars_free(Sol1);
#endif
    
    ElemSpaceDiscr_free(elem);
    NodeSpaceDiscr_free(node);
}




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: data.c,v $
 Revision 1.2  1998/03/07 09:56:45  nicola
 Added flux computation and multiple pressure variables.
 
 Revision 1.1  1998/03/01 18:43:33  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
