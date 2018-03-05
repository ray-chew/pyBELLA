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
 Mupltiple pressure variables embedded in a Multigrid structure
 ------------------------------------------------------------------------------*/
typedef struct {
	
    int level_number;
	
    ElemSpaceDiscr *elem;
    NodeSpaceDiscr *node;
    
    double *p_coarse;
    double *p;
    double *rhs;
    double *diaginv;
	
    double *wplus[3];
    /*
    double *wminus[3];
    double *lwplus[3];
    double *lwminus[3];
     */
	double *wcenter;
	double *wgrav;
	
} LowMachMGLevel;


typedef struct {
	
	double p0;
	double p00;
    double *p2_cells;    
    double *dp2_cells;    
    double *p2_nodes;
    double *dp2_nodes;
    double *eta_hyp0;
    double *eta_hyp;
    double *eta;
    double *eta0;
	
    States *HydroState;
    States *HydroState_n;
	
	double dt;
    LowMachMGLevel** Level;
	
	/* correction velocities for open boundaries */
	double du;
	double Sdu;
	
} MPV;



/*------------------------------------------------------------------------------
 allocation / freeing
 ------------------------------------------------------------------------------*/
void initialize_projection(
							  const int inx,
							  const int iny,
							  const int inz,
							  const double x0, 
							  const double x1, 
							  const double y0, 
							  const double y1, 
							  const double z0, 
							  const double z1,
							  const enum BdryType left, 
							  const enum BdryType right, 
							  const enum BdryType bottom, 
							  const enum BdryType top, 
							  const enum BdryType back, 
							  const enum BdryType front);

void terminate_MG_projection(
							 void);

void initialize_HydroState(
						   int nc,
                           int nn);

MPV* MPV_new(const Grid* base_grid);

void MPV_free(
			  MPV* var);

LowMachMGLevel* LowMachMGLevel_new(
								   const int level, 
								   const Grid* base_grid);

void LowMachMGLevel_free(
						 LowMachMGLevel* var);

#endif /* MPV_H */

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: mpv.h,v $
 Revision 1.1  1998/03/07 09:56:47  nicola
 Added flux computation and multiple pressure variables.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
