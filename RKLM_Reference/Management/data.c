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
#include "float.h"
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

/* Context */
int ndim;
int H1;
int H2;
int H3; 

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
ConsVars* dSol;           /* full size */ /* TODO: Can I work without full-size dSol arrays? */
VectorField* adv_flux;    /* full size, components located on primary cell faces */
ConsVars* flux[3];        /* full size (M < 1.0) */

#ifdef SYMMETRIC_ADVECTION
ConsVars* Sol1;           /* full size (M < 1.0; needed to test x-y-symmetric advection for the travelling vortex) */
#endif

States* Solk;             /* cache size */
double* W0;               /* full nG-length double array */


void Data_init() {
	
	int i, inx, iny, inz, n_aux;
	double x0, x1, y0, y1, z0, z1;
	enum BdryType left, right, bottom, top, back, front;
	Grid* grid;
	
	/* User data */
	User_Data_init(&ud);
	
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
	
	/* Context */
	ndim = inz > 1 ? 3 : (iny > 1 ? 2 : 1);
	H1 = HEAVISIDE(ndim);
	H2 = HEAVISIDE(ndim - 1);
	H3 = HEAVISIDE(ndim - 2);
	
	/* Time discretization */
	step = 0;
	t = 0.0;
	
	/* Grid and space discretization */
	grid = Grid_new(
					inx, 
					iny, 
					inz, 
					x0, 
					x1, 
					y0, 
					y1, 
					z0, 
					z1, 
					left,
					right,
					bottom,
					top,
					back,
					front); 
	elem = ElemSpaceDiscr_new(grid);
	node = NodeSpaceDiscr_new(grid);
	Grid_free(grid); 
	
	/* Thermodynamics */
	Thermodynamic_init(&th, &ud);
	
	/* Arrays */
	Sol  = ConsVars_new(elem->nc);
	dSol = ConsVars_new(elem->nc);
    adv_flux  = VectorField_new(node->nc);
	Solk = States_small_new(3 * ud.ncache / 2);
    
    /* auxiliary space */
    n_aux = node->ifx;
    n_aux *= (node->ndim > 1 ? node->ify : 1);
    n_aux *= (node->ndim > 2 ? node->ifz : 1);
    
    W0  = (double*)malloc((unsigned)(n_aux * sizeof(double)));
    
    Sol0    = ConsVars_new(elem->nc);
#ifdef SYMMETRIC_ADVECTION
    Sol1    = ConsVars_new(elem->nc);
#endif
    
    flux[0] = ConsVars_new(elem->nfx);
    ConsVars_setzero(flux[0], elem->nfx);
    if(ndim > 1) {
        flux[1] = ConsVars_new(elem->nfy);
        ConsVars_setzero(flux[1], elem->nfy);
    }
    if(ndim > 2) {
        flux[2] = ConsVars_new(elem->nfz);
        ConsVars_setzero(flux[2], elem->nfz);
    }
    
    initialize_bdry(elem);
    
    initialize_projection(inx,iny,inz,
                          x0, x1, y0, y1, z0, z1,
                          left,right,bottom,top,back,front);
    
    
    for (i=0; i<3; i++) {
		if (ud.i_gravity[i]) initialize_HydroState(elem->ic[i], node->ic[i]);
	}
	
	Sol_initial(Sol, elem, node);
	ConsVars_set(Sol0, Sol, elem->nc);
	ConsVars_setzero(dSol, elem->nc);

#ifdef SYMMETRIC_ADVECTION
    ConsVars_setzero(Sol1, elem->nc);
#endif
    
    Explicit_malloc(3 * ud.ncache / 2);
	recovery_malloc(3 * ud.ncache / 2);  

}


void Data_free() {
	
	ElemSpaceDiscr_free(elem);
	NodeSpaceDiscr_free(node);
	ConsVars_free(Sol);
	ConsVars_free(dSol);
    VectorField_free(adv_flux);
	States_small_free(Solk);
	free(W0);
    ConsVars_free(Sol0);

    if(ndim > 2) ConsVars_free(flux[2]);
    if(ndim > 1) ConsVars_free(flux[1]);
    ConsVars_free(flux[0]);

    terminate_MG_projection();
    
    Explicit_free();
    recovery_free();    
    close_bdry();
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
