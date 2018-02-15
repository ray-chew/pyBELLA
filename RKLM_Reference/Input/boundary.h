/*******************************************************************************
 File:   boundary.h
 Author: Nicola 
 Date:   Thu Feb 19 07:16:28 CET 1998
 *******************************************************************************/
#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "enumerator.h"
#include "variable.h"
#include "mpv.h"

typedef struct {
	
	double* wall_massflux;
	double* wall_slope;
	double* wall_relative_slope;
} BDRY;

/*------------------------------------------------------------------------------
 
 ------------------------------------------------------------------------------*/

void initialize_bdry(
					 const ElemSpaceDiscr* elem);

void close_bdry( 
                void );

void Bound(
		   ConsVars* Sol, 
		   const double lambda, 
		   const int n, 
           const int SplitStep, 
           const int setZ);

void set_wall_massflux(
					   BDRY* bdry, 
					   const ConsVars* Sol0, 
					   const ElemSpaceDiscr* elem);

double slanted_wall_slope(
						  double x);

double velo_background(
					   double t);

double wall_massflux(
					 double x,
					 double y,
					 double wind_speed_x,
					 double wind_speed_y);

void check_flux_bcs(
					States* Lefts, 
					States* Rights,
					const int nmax,
					const int kcache,
					const int njump,
					ElemSpaceDiscr* elem,
					const int SplitStep);


/* ============================================================================= */

void Set_Explicit_Boundary_Data(
                                ConsVars* Sol,
                                const ElemSpaceDiscr* elem,
                                const MPV* mpv,
                                const int setZ);


/*
 void check_bdry_fluxes(
 ConsVars* Fluxes,
 const int nmax,
 const int kcache,
 const int njump,
 ElemSpaceDiscr* elem,
 const int SplitStep);
 */				    
#endif /* BOUNDARY_H */




/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: boundary.h,v $
 Revision 1.1  1998/03/01 18:43:33  nicola
 This is the initial revision of 3d. It comes out after two weeks of work on
 Matthias' version of Rupert's F3D code. It is the starting point for imple
 menting a low Mach number extension.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
