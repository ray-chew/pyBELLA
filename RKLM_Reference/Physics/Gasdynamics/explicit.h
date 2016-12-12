/*******************************************************************************
 File:   explicit.h
 Author: Nicola, Rupert
 Date:   Mon Mar  2 10:47:26 CET 1998 / 2013.03.11
 *******************************************************************************/
#ifndef EXPLICIT_H
#define EXPLICIT_H

#include "variable.h"
#include "mpv.h"

/*------------------------------------------------------------------------------
 Parameters for the predictor step time integrator
 ------------------------------------------------------------------------------*/
#define NO_OF_RK_STAGES 3
typedef struct {
    double dt_frac;
    double flux_frac[NO_OF_RK_STAGES][2];
    double update_frac[NO_OF_RK_STAGES];
    int no_of_stages;
    enum Boolean multiD_updt;
} TimeIntegratorParams;


/*------------------------------------------------------------------------------
 Allocates memory for the arrays used by Explicit_step() and  Explicit_flux()
 ------------------------------------------------------------------------------*/
void Explicit_malloc(const int size);


/*------------------------------------------------------------------------------
 Releases memory for the arrays used by Explicit_step() and  Explicit_flux()
 ------------------------------------------------------------------------------*/
void Explicit_free();


/*------------------------------------------------------------------------------
 Explicit step and flux computation
 ------------------------------------------------------------------------------*/
void Explicit_step_and_flux(
							ConsVars* Sol,
							ConsVars* flux,
                            double* buoyS,
                            VectorField* buoy,
                            double* dp2,
							const States* HydroState,
							const double lambda, 
							const int n, 
							const int SplitStep,
                            const int RK_stage,
                            const int implicit);


/*------------------------------------------------------------------------------
 Pressure and gravity in split mode
 ------------------------------------------------------------------------------*/
void Explicit_pressure_and_gravity(
                                   ConsVars* Sol,
                                   double* buoyS,
                                   VectorField* buoy,
                                   const MPV* mpv,
                                   const ElemSpaceDiscr* elem,
                                   const NodeSpaceDiscr* node,
                                   const double dt);

/*------------------------------------------------------------------------------
 Boundary Absorber
 ------------------------------------------------------------------------------*/
void Absorber(
			  ConsVars* Sol,
			  double time,
			  double dt);

/*------------------------------------------------------------------------------
 Update
 ------------------------------------------------------------------------------*/
void Explicit_step_update( 
						  ConsVars* Sol,
						  const int n);



/*------------------------------------------------------------------------------
 Computes updates from full-dimensional flux balances
 ------------------------------------------------------------------------------*/
void fullD_explicit_updates(ConsVars* Sol, 
                            ConsVars* Sol0,
                            ConsVars* flux[3],
                            double* buoyS,
                            VectorField* buoy, 
                            const ElemSpaceDiscr* elem, 
                            const double dt,
                            const int RK_stage);

void Explicit_Coriolis(ConsVars *Sol, const ElemSpaceDiscr* elem, const double dt);


#ifdef GRAVITY_IMPLICIT_2
/*------------------------------------------------------------------------------
 explicit step for the fast linear system
 ------------------------------------------------------------------------------*/
void Explicit_Buoyancy(ConsVars* Sol, 
                       VectorField* buoy, 
                       const MPV* mpv, 
                       const ElemSpaceDiscr* elem, 
                       const NodeSpaceDiscr* node, 
                       const double t, 
                       const double dt,
                       const int implicit);
#endif

#endif /* EXPLICIT_H */

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: explicit.h,v $
 Revision 1.1  1998/03/07 09:56:46  nicola
 Added flux computation and multiple pressure variables.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
