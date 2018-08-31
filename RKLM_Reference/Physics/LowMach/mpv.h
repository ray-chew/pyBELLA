/*******************************************************************************
 File:   mpv.h
 Author: Nicola
 Date:   Fri Mar  6 16:06:19 WET 1998
 *******************************************************************************/
#ifndef MPV_H
#define MPV_H

#include <assert.h> 

#include "kgrid.h"
#include "variable.h"
/* #include "space_discretization.h" */

/*------------------------------------------------------------------------------
 Structure needed for semi-implicit time integration
 ------------------------------------------------------------------------------*/
typedef struct {
	
	double p0;
	double p00;
    double *p2_cells;    
    double *dp2_cells;    
    double *p2_nodes;
    double *dp2_nodes;

    double *drhou_cells;
    double *drhov_cells;
    double *drhow_cells;

    double *rhs;
    double *diaginv;
    
    double *wplus[3];
    double *wcenter;
    double *wgrav;
    
    States *HydroState;
    States *HydroState_n;
	
	double dt;

	/* correction velocities for open boundaries */
	double du;
	double Sdu;
	
} MPV;

/*------------------------------------------------------------------------------
 allocation / freeing
 ------------------------------------------------------------------------------*/
MPV* MPV_new(
             const ElemSpaceDiscr *elem,
             const NodeSpaceDiscr *node);

void MPV_free(
			  MPV* var,
              const ElemSpaceDiscr *elem);

#endif /* MPV_H */

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: mpv.h,v $
 Revision 1.1  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
