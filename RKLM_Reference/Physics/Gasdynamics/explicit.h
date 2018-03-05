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
    enum Boolean multiD_updt;
} TimeIntegratorParams;


/*------------------------------------------------------------------------------
 Allocates memory for the arrays used by Explicit_step() and  Explicit_flux()
 ------------------------------------------------------------------------------*/
void Explicit_malloc(const int size);


/*------------------------------------------------------------------------------
 Releases memory for the arrays used by Explicit_step() and  Explicit_flux()
 ------------------------------------------------------------------------------*/
void Explicit_free( void );


/*------------------------------------------------------------------------------
 Advection cast into a single call routine
 ------------------------------------------------------------------------------*/
void advect(
            ConsVars *Sol, 
            ConsVars* flux[3],
            const double dt, 
            const ElemSpaceDiscr* elem,
            const enum FluxesFrom adv_fluxes_from, 
            const enum MUSCL_ON_OFF muscl_on_off, 
            const enum No_of_Strang_Sweeps no_of_sweeps);

/*------------------------------------------------------------------------------
 Explicit step and flux computation
 ------------------------------------------------------------------------------*/
void Explicit_step_and_flux(
                            ConsVars* Sol,
                            ConsVars* flux,
                            const double lambda, 
                            const int n, 
                            const int SplitStep,
                            const int RK_stage,
                            const enum FluxesFrom adv_fluxes_from, 
                            const enum MUSCL_ON_OFF muscl_on_off);


/*------------------------------------------------------------------------------
 Pressure and gravity in split mode
 ------------------------------------------------------------------------------*/
void Explicit_pressure_and_gravity(
                                   ConsVars* Sol,
                                   const MPV* mpv,
                                   const ElemSpaceDiscr* elem,
                                   const NodeSpaceDiscr* node,
                                   const double dt);

/*------------------------------------------------------------------------------
 Boundary Absorber
 ------------------------------------------------------------------------------*/
void Absorber(
              ConsVars* Sol,
              const ElemSpaceDiscr* elem,
              const double time,
              const double dt);

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
                            const ElemSpaceDiscr* elem, 
                            const double dt,
                            const int RK_stage) ;

void Explicit_Coriolis(ConsVars *Sol, const ElemSpaceDiscr* elem, const double dt);

#endif /* EXPLICIT_H */

/*LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
 $Log: explicit.h,v $
 Revision 1.1  1998/03/07 09:56:46  nicola
 Added flux computation and multiple pressure variables.
 
 LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL*/
